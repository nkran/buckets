"""FALLBACK generator: build valid Escher map JSON directly from iCHO3K.

The primary maps come from curated RECON1 human pathway maps (see
`scripts/prepare_escher_maps.py`), which are real, connected, hand-laid-out cascades
and overlay iCHO3K at 73–88% reaction coverage. This generator is a fallback for
subsystems that lack a curated map: it guarantees 100% id match (maps are built from
the model) but only produces a functional "reaction-card" grid — each reaction gets a
midmarker with its (non-currency) reactants on the left and products on the right,
NOT a connected pathway. Use curated maps whenever available.

Escher map schema (v1-0-0): a two-element array
    [ {schema, homepage, map_id, map_name, map_description},
      {nodes, reactions, text_labels, canvas} ]
"""

from __future__ import annotations

import re
from pathlib import Path

SCHEMA = "https://escher.github.io/escher/jsonschema/1-0-0#"
HOMEPAGE = "https://escher.github.io"

# currency / cofactor metabolites hidden from the layout to reduce clutter
CURRENCY = {
    "h", "h2o", "atp", "adp", "amp", "pi", "ppi", "nad", "nadh", "nadp", "nadph",
    "co2", "o2", "coa", "nh4", "hco3", "fad", "fadh2", "fmn", "gtp", "gdp", "gmp",
    "so4", "h2o2", "na1", "k", "cl", "ca2", "mg2", "fe2", "fe3", "o2s", "udp", "ump",
    "utp", "ctp", "cdp", "cmp", "accoa", "q10", "q10h2", "pmtcoa",
}

_COMPARTMENT_RE = re.compile(r"_([a-z]{1,2})$")


def base_metabolite_id(mid: str) -> str:
    """Strip a trailing compartment suffix (e.g. `atp_c` -> `atp`)."""
    return _COMPARTMENT_RE.sub("", mid)


def _bezier(p_from: tuple[float, float], p_to: tuple[float, float]):
    """Two control points at 1/3 and 2/3 along the segment (straight-ish arrows)."""
    (x0, y0), (x1, y1) = p_from, p_to
    return ({"x": x0 + (x1 - x0) / 3, "y": y0 + (y1 - y0) / 3},
            {"x": x0 + 2 * (x1 - x0) / 3, "y": y0 + 2 * (y1 - y0) / 3})


def generate_subsystem_map(model, subsystem: str, *, subsystems: set[str] | None = None,
                           max_reactions: int = 60, cols: int = 5, dx: float = 420,
                           dy: float = 320, met_gap: float = 62,
                           hide_currency: bool = True) -> list:
    """Build an Escher map (as a Python list, ready for json.dump).

    Selects reactions whose subsystem == ``subsystem``, or whose subsystem is in
    ``subsystems`` (a set) when given — in which case ``subsystem`` is used only as
    the map title (e.g. a "Central carbon" overview spanning several subsystems).
    """
    wanted = subsystems if subsystems is not None else {subsystem}
    rxns = [r for r in model.reactions if (r.subsystem or "").strip() in wanted]
    truncated = len(rxns) > max_reactions
    rxns = rxns[:max_reactions]

    nodes: dict[str, dict] = {}
    reactions: dict[str, dict] = {}
    text_labels: dict[str, dict] = {}
    xs: list[float] = []
    ys: list[float] = []
    counter = {"n": 0}

    def new_id() -> str:
        counter["n"] += 1
        return str(counter["n"])

    def add_metabolite_node(met, x, y, primary) -> str:
        nid = new_id()
        nodes[nid] = {"node_type": "metabolite", "x": x, "y": y,
                      "bigg_id": met.id, "name": met.name or met.id,
                      "label_x": x + 10, "label_y": y - 10,
                      "node_is_primary": bool(primary)}
        xs.append(x); ys.append(y)
        return nid

    for i, r in enumerate(rxns):
        col, row = i % cols, i // cols
        cx, cy = col * dx + 220, row * dy + 160

        mid_id = new_id()
        nodes[mid_id] = {"node_type": "midmarker", "x": cx, "y": cy}
        xs.append(cx); ys.append(cy)

        mets = list(r.metabolites.items())
        disp = mets
        if hide_currency:
            filtered = [(m, c) for m, c in mets if base_metabolite_id(m.id) not in CURRENCY]
            disp = filtered or mets
        reactants = [(m, c) for m, c in disp if c < 0][:4]
        products = [(m, c) for m, c in disp if c > 0][:4]

        segments: dict[str, dict] = {}
        used_mets: list[tuple] = []

        def place_side(items, sign):
            n = len(items)
            for k, (m, c) in enumerate(items):
                mx = cx + sign * 150
                my = cy + (k - (n - 1) / 2) * met_gap
                node = add_metabolite_node(m, mx, my, primary=(n == 1))
                if sign < 0:      # reactant -> midmarker
                    b1, b2 = _bezier((mx, my), (cx, cy))
                    segments[new_id()] = {"from_node_id": node, "to_node_id": mid_id,
                                          "b1": b1, "b2": b2}
                else:             # midmarker -> product
                    b1, b2 = _bezier((cx, cy), (mx, my))
                    segments[new_id()] = {"from_node_id": mid_id, "to_node_id": node,
                                          "b1": b1, "b2": b2}
                used_mets.append({"coefficient": float(c), "bigg_id": m.id})

        place_side(reactants, -1)
        place_side(products, +1)

        reactions[new_id()] = {
            "name": r.name or r.id,
            "bigg_id": r.id,
            "reversibility": bool(r.reversibility),
            "label_x": cx + 10, "label_y": cy - 25,
            "gene_reaction_rule": r.gene_reaction_rule or "",
            "genes": [{"bigg_id": g.id, "name": g.name or g.id} for g in r.genes],
            "metabolites": used_mets,
            "segments": segments,
        }

    if not xs:  # empty subsystem — produce a minimal valid empty map
        xs, ys = [0.0], [0.0]
    title = new_id()
    text_labels[title] = {"x": min(xs), "y": min(ys) - 60,
                          "text": subsystem + (" (truncated)" if truncated else "")}

    pad = 120
    canvas = {"x": min(xs) - pad, "y": min(ys) - pad,
              "width": (max(xs) - min(xs)) + 2 * pad,
              "height": (max(ys) - min(ys)) + 2 * pad}

    meta = {"schema": SCHEMA, "homepage": HOMEPAGE,
            "map_id": re.sub(r"[^A-Za-z0-9]+", "_", subsystem).strip("_"),
            "map_name": subsystem,
            "map_description": f"iCHO3K {subsystem} — auto-generated"
                               + (" (reactions truncated)" if truncated else "")}
    return [meta, {"nodes": nodes, "reactions": reactions,
                   "text_labels": text_labels, "canvas": canvas}]


def map_reaction_ids(escher_map: list) -> set[str]:
    """Return the set of reaction BiGG ids present in an Escher map."""
    return {r["bigg_id"] for r in escher_map[1]["reactions"].values()}


def write_map(escher_map: list, path: Path) -> Path:
    import json
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(escher_map))
    return path
