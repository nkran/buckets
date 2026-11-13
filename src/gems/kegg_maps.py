"""Assemble connected Escher pathway maps from multiple sources.

Layout + connectivity come from **human (hsa) KEGG KGML** (hand-drawn coordinates +
substrate/product links, restricted to reactions that occur in humans — a good proxy
for what the Recon3D-derived iCHO3K should contain). Reaction identities are bridged to
iCHO3K through a **union crosswalk**: the BiGG KEGG↔BiGG links *plus* a MetaNetX (MNX)
bridge (``reac_xref``), which recovers reactions BiGG doesn't link to KEGG directly.

Readability improvements over raw KGML: metabolites and reactions carry human-readable
names, and currency metabolites (ATP, H2O, NAD, …) are drawn as secondary (small) nodes
so the pathway backbone stands out. ``coverage`` (mapped/total pathway reactions) is the
headline reliability metric.
"""

from __future__ import annotations

import re
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from . import config as C
from . import escher_maps as em

BIGG_RXN_NS = C.DATA / "maps" / "bigg_models_reactions.txt"
BIGG_MET_NS = C.DATA / "maps" / "bigg_models_metabolites.txt"
MNX_XREF = C.DATA / "maps" / "mnx_reac_xref.tsv"
KEGG_CPD_NAMES = C.DATA / "maps" / "kegg_compound_names.txt"
KGML_URL = "https://rest.kegg.jp/get/{pid}/kgml"

_COMPARTMENTS = "mcxrnglei"

# KEGG currency/cofactor compounds -> drawn as secondary (small) nodes
CURRENCY_KEGG = {
    "C00001", "C00002", "C00003", "C00004", "C00005", "C00006", "C00007", "C00008",
    "C00009", "C00010", "C00011", "C00013", "C00014", "C00015", "C00016", "C00019",
    "C00020", "C00021", "C00035", "C00044", "C00055", "C00075", "C00080", "C00131",
    "C00144", "C00360", "C00390", "C01352", "C00028", "C00030", "C00006", "C00342",
}


def _strip_compartment(bigg_id: str) -> str:
    return re.sub(rf"[{_COMPARTMENTS}]$", "", bigg_id)


# ---------------------------------------------------------------------------
# Crosswalks
# ---------------------------------------------------------------------------
def load_reaction_crosswalk(model) -> tuple[dict[str, set[str]], dict[str, int]]:
    """KEGG reaction Rxxxxx -> set(iCHO3K reaction ids), unioning two bridges:

    1. direct BiGG ``kegg.reaction`` links, and
    2. a MetaNetX bridge from ``reac_xref`` (KEGG->MNX and BiGG->MNX in one version).

    Returns (kegg2model, source_counts) where source_counts summarises recall.
    """
    model_rxn = {r.id for r in model.reactions}
    key_index: dict[str, set[str]] = defaultdict(set)
    for rid in model_rxn:
        key_index[rid].add(rid)
        key_index[_strip_compartment(rid)].add(rid)

    kegg2model: dict[str, set[str]] = defaultdict(set)

    # (1) direct BiGG kegg.reaction links
    n_direct = 0
    with open(BIGG_RXN_NS) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            krs = re.findall(r"kegg\.reaction/(R\d+)", p[4])
            if not krs:
                continue
            hits = key_index.get(p[0], set()) | key_index.get(_strip_compartment(p[0]), set())
            for kr in krs:
                if hits:
                    kegg2model[kr] |= hits
                    n_direct += 1

    # (2) MetaNetX bridge (reac_xref: kegg.reaction:R -> MNX ; bigg.reaction:ID -> MNX)
    n_mnx = 0
    if MNX_XREF.exists():
        mnx2kegg: dict[str, set[str]] = defaultdict(set)
        bigg2mnx: dict[str, set[str]] = defaultdict(set)
        with open(MNX_XREF) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2 or not parts[1].startswith("MNXR"):
                    continue
                src, mnx = parts[0], parts[1]
                if src.startswith("kegg.reaction:"):
                    mnx2kegg[mnx].add(src.split(":", 1)[1])
                elif src.startswith("bigg.reaction:"):
                    bigg2mnx[src.split(":", 1)[1]].add(mnx)
        mnx2model: dict[str, set[str]] = defaultdict(set)
        for rid in model_rxn:
            for m in (bigg2mnx.get(rid) or bigg2mnx.get(_strip_compartment(rid)) or ()):
                mnx2model[m].add(rid)
        for mnx, krs in mnx2kegg.items():
            mm = mnx2model.get(mnx)
            if not mm:
                continue
            for kr in krs:
                if kr not in kegg2model or not (kegg2model[kr] & mm):
                    n_mnx += 1
                kegg2model[kr] |= mm

    kegg2model = {k: v for k, v in kegg2model.items() if v}
    return kegg2model, {"direct": n_direct, "mnx_added": n_mnx, "total": len(kegg2model)}


def load_compound_names() -> dict[str, str]:
    names: dict[str, str] = {}
    if KEGG_CPD_NAMES.exists():
        for line in KEGG_CPD_NAMES.read_text().splitlines():
            cid, _, rest = line.partition("\t")
            cid = cid.replace("cpd:", "")
            if rest:
                names[cid] = rest.split(";")[0].strip()
    return names


def load_metabolite_crosswalk() -> dict[str, str]:
    """KEGG compound Cxxxxx -> a BiGG metabolite base id (for overlay/labels)."""
    c2bigg: dict[str, str] = {}
    if BIGG_MET_NS.exists():
        with open(BIGG_MET_NS) as fh:
            next(fh)
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                for c in re.findall(r"kegg\.compound/(C\d+)", p[4]):
                    c2bigg.setdefault(c, p[0])
    return c2bigg


def build_reaction_names(model) -> dict[str, str]:
    return {r.id: (r.name or r.id) for r in model.reactions}


# ---------------------------------------------------------------------------
# KGML parsing + map building
# ---------------------------------------------------------------------------
@dataclass
class MapStats:
    pathway_id: str
    title: str
    n_reactions_total: int = 0
    n_reactions_mapped: int = 0
    n_ambiguous: int = 0
    # displayed reaction id -> all iCHO3K reactions matching that KEGG step (when >1)
    aliases: dict = field(default_factory=dict)

    @property
    def coverage(self) -> float:
        return self.n_reactions_mapped / max(self.n_reactions_total, 1)


def fetch_kgml(pathway_id: str, cache_dir: Path = C.DATA / "maps" / "kgml") -> str:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = cache_dir / f"{pathway_id}.xml"
    if cached.exists():
        return cached.read_text()
    req = urllib.request.Request(KGML_URL.format(pid=pathway_id),
                                 headers={"User-Agent": "gems/0.1"})
    text = urllib.request.urlopen(req, timeout=30).read().decode()
    cached.write_text(text)
    return text


def build_map_from_kgml(pathway_id: str, kgml_text: str, kegg2model: dict[str, set[str]],
                        c2bigg: dict[str, str], cpd_names: dict[str, str],
                        rxn_names: dict[str, str], *, scale: float = 1.7) -> tuple[list, MapStats]:
    root = ET.fromstring(kgml_text)
    title = re.sub(r"\s*-\s*Homo sapiens.*$", "", root.get("title", pathway_id))

    comp_entries: dict[str, tuple[str, float, float]] = {}
    for e in root.findall("entry"):
        if e.get("type") != "compound":
            continue
        g = e.find("graphics")
        if g is None or g.get("x") is None:
            continue
        kegg_c = e.get("name", "").split()[0].replace("cpd:", "")
        comp_entries[e.get("id")] = (kegg_c, float(g.get("x")), float(g.get("y")))

    rxn_coords: dict[str, tuple[float, float]] = {}
    for e in root.findall("entry"):
        if e.get("type") not in ("reaction", "gene"):
            continue
        g = e.find("graphics")
        if g is None or g.get("x") is None:
            continue
        for nm in (e.get("reaction") or "").split():
            if nm.startswith("rn:"):
                rxn_coords.setdefault(nm[3:], (float(g.get("x")), float(g.get("y"))))

    nodes: dict[str, dict] = {}
    reactions: dict[str, dict] = {}
    xs: list[float] = []
    ys: list[float] = []
    counter = {"n": 0}
    node_for_entry: dict[str, str] = {}

    def nid() -> str:
        counter["n"] += 1
        return str(counter["n"])

    def add_compound_node(entry_id: str) -> str | None:
        if entry_id in node_for_entry:
            return node_for_entry[entry_id]
        if entry_id not in comp_entries:
            return None
        kegg_c, x, y = comp_entries[entry_id]
        x, y = x * scale, y * scale
        bigg = c2bigg.get(kegg_c, kegg_c)
        name = cpd_names.get(kegg_c, bigg)
        primary = kegg_c not in CURRENCY_KEGG
        n = nid()
        nodes[n] = {"node_type": "metabolite", "x": x, "y": y, "bigg_id": bigg,
                    "name": name, "label_x": x + 12, "label_y": y + 10,
                    "node_is_primary": primary}
        xs.append(x); ys.append(y)
        node_for_entry[entry_id] = n
        return n

    stats = MapStats(pathway_id=pathway_id, title=title)
    seen: set[str] = set()

    for rel in root.findall("reaction"):
        kegg_rs = [nm[3:] for nm in rel.get("name", "").split() if nm.startswith("rn:")]
        mapped_ids: set[str] = set()
        for kr in kegg_rs:
            mapped_ids |= kegg2model.get(kr, set())
        for kr in kegg_rs:
            if kr in seen:
                continue
            seen.add(kr)
            stats.n_reactions_total += 1
            hits = kegg2model.get(kr, set())
            if hits:
                stats.n_reactions_mapped += 1
            if len(hits) > 1:
                stats.n_ambiguous += 1
        if not mapped_ids:
            continue
        bigg_rxn = sorted(mapped_ids)[0]
        if len(mapped_ids) > 1:
            stats.aliases[bigg_rxn] = sorted(mapped_ids)

        coord = next((rxn_coords[kr] for kr in kegg_rs if kr in rxn_coords), None)
        subs = [s.get("id") for s in rel.findall("substrate")]
        prods = [p.get("id") for p in rel.findall("product")]
        sub_nodes = [n for n in (add_compound_node(s) for s in subs) if n]
        prod_nodes = [n for n in (add_compound_node(p) for p in prods) if n]
        if not sub_nodes and not prod_nodes:
            continue
        if coord is None:
            pts = [(nodes[n]["x"], nodes[n]["y"]) for n in sub_nodes + prod_nodes]
            coord = (sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts))
        else:
            coord = (coord[0] * scale, coord[1] * scale)
        mid = nid()
        nodes[mid] = {"node_type": "midmarker", "x": coord[0], "y": coord[1]}
        xs.append(coord[0]); ys.append(coord[1])

        segments: dict[str, dict] = {}
        for n in sub_nodes:
            b1, b2 = em._bezier((nodes[n]["x"], nodes[n]["y"]), coord)
            segments[nid()] = {"from_node_id": n, "to_node_id": mid, "b1": b1, "b2": b2}
        for n in prod_nodes:
            b1, b2 = em._bezier(coord, (nodes[n]["x"], nodes[n]["y"]))
            segments[nid()] = {"from_node_id": mid, "to_node_id": n, "b1": b1, "b2": b2}

        reactions[nid()] = {
            "name": rxn_names.get(bigg_rxn, bigg_rxn), "bigg_id": bigg_rxn,
            "reversibility": rel.get("type") == "reversible",
            "label_x": coord[0] + 12, "label_y": coord[1] - 12,
            "gene_reaction_rule": "", "genes": [],
            "metabolites": ([{"coefficient": -1.0, "bigg_id": nodes[n]["bigg_id"]} for n in sub_nodes]
                            + [{"coefficient": 1.0, "bigg_id": nodes[n]["bigg_id"]} for n in prod_nodes]),
            "segments": segments,
        }

    if not xs:
        xs, ys = [0.0], [0.0]
    pad = 100
    canvas = {"x": min(xs) - pad, "y": min(ys) - pad,
              "width": (max(xs) - min(xs)) + 2 * pad, "height": (max(ys) - min(ys)) + 2 * pad}
    meta = {"schema": em.SCHEMA, "homepage": em.HOMEPAGE, "map_id": pathway_id,
            "map_name": title,
            "map_description": f"KEGG {pathway_id} (human) → iCHO3K "
                               f"({stats.n_reactions_mapped}/{stats.n_reactions_total} reactions)"}
    return [meta, {"nodes": nodes, "reactions": reactions, "text_labels": {}, "canvas": canvas}], stats
