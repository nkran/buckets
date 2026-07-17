"""Assemble as many Escher pathway maps as possible, from multiple sources, and
score each one's reliability so overlays can be interpreted correctly.

Sources
-------
1. **Curated RECON1 human maps** (escher.github.io) — hand-drawn, connected, use BiGG
   ids directly; overlay iCHO3K at 73–88%. Highest reliability.
2. **KEGG-derived maps** — layout + connectivity from KEGG KGML, reaction identities
   bridged to iCHO3K through the BiGG KEGG↔BiGG crosswalk. Covers many more pathways;
   reliability varies with how well each pathway maps (reported per map).

Reliability tiers (written into `escher_maps/manifest.json`):
  * **High**   – curated RECON1 map (direct BiGG ids).
  * **Medium** – KEGG map, ≥40% of pathway reactions mapped to iCHO3K.
  * **Low**    – KEGG map, <40% mapped (core reactions shown, many gaps).
Each entry also carries `coverage`, `n_mapped`, `n_reactions`, and `ambiguous_frac`
(fraction of mapped reactions where one KEGG reaction matched several iCHO3K reactions
and one was chosen).

Run:  python scripts/prepare_escher_maps.py
"""

from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

from gems import config as C
from gems import kegg_maps as km
from gems import model_setup

MAPS_DIR = C.ROOT / "escher_maps"
ESCHER_JS_URL = "https://unpkg.com/escher@1.8.1/dist/escher.min.js"
ESCHER_MAP_BASE = "https://escher.github.io/1-0-0/6/maps"
KEGG_LIST_URL = "https://rest.kegg.jp/list/pathway/hsa"   # human-specific pathways

MIN_MAPPED = 5           # skip KEGG maps with fewer than this many iCHO3K reactions
MEDIUM_COVERAGE = 0.40   # tier threshold

RECON1 = [
    ("Homo sapiens", "RECON1.Glycolysis TCA PPP", "Central carbon", "glycolysis_tca_ppp"),
    ("Homo sapiens", "RECON1.Carbohydrate metabolism", "Carbohydrate", "carbohydrate_metabolism"),
    ("Homo sapiens", "RECON1.Amino acid metabolism (partial)", "Amino acid", "amino_acid_metabolism"),
    ("Homo sapiens", "RECON1.Tryptophan metabolism", "Amino acid", "tryptophan_metabolism"),
    ("Homo sapiens", "RECON1.Inositol retinol metabolism", "Lipid & vitamins", "inositol_retinol_metabolism"),
]


def _download(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "gems/0.1"})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return fh.read()


def vendor_escher_js() -> None:
    dest = MAPS_DIR / "lib" / "escher.min.js"
    if dest.exists() and dest.stat().st_size > 100_000:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        dest.write_bytes(_download(ESCHER_JS_URL))
        print(f"vendored escher.min.js ({dest.stat().st_size} bytes)")
    except Exception as exc:
        print(f"WARNING: could not download escher.min.js ({exc})")


def kegg_category(num: int) -> str:
    """Broad category from the KEGG pathway number (standard KEGG ranges)."""
    def inr(*rngs):
        return any(a <= num <= b for a, b in rngs)
    if inr((10, 53), (500, 562), (620, 660)):
        return "Carbohydrate"
    if inr((190, 199), (680, 680), (710, 720), (910, 920)):
        return "Energy"
    if inr((61, 73), (100, 140), (561, 565), (590, 592), (600, 600), (1040, 1040)):
        return "Lipid"
    if inr((230, 240)):
        return "Nucleotide"
    if inr((220, 220), (250, 300), (310, 360), (380, 380), (400, 400), (470, 480)):
        return "Amino acid"
    if inr((510, 515), (520, 520), (531, 534), (563, 563), (601, 604)):
        return "Glycan"
    if inr((130, 130), (670, 670), (730, 790), (830, 830), (860, 860), (900, 900)):
        return "Cofactors & vitamins"
    return "Other"


def fetch_recon1(model_rxn: set[str], manifest: dict) -> None:
    for organism, map_name, category, slug in RECON1:
        url = f"{ESCHER_MAP_BASE}/{urllib.parse.quote(organism)}/{urllib.parse.quote(map_name)}.json"
        try:
            emap = json.loads(_download(url))
        except Exception as exc:
            print(f"SKIP RECON1 {map_name}: {exc}")
            continue
        map_rxn = {r["bigg_id"] for r in emap[1]["reactions"].values()}
        overlap = map_rxn & model_rxn
        coverage = len(overlap) / max(len(map_rxn), 1)
        out = MAPS_DIR / f"{slug}.json"
        out.write_text(json.dumps(emap))
        manifest.setdefault(category, []).append({
            "pathway": emap[0].get("map_name", map_name), "file": out.name,
            "source": "RECON1 (curated)", "tier": "High",
            "coverage": round(coverage, 3), "n_reactions": len(map_rxn),
            "n_mapped": len(overlap), "ambiguous_frac": 0.0,
        })
        print(f"[RECON1] {map_name:<40} {len(overlap)}/{len(map_rxn)} ({coverage:.0%})")


def fetch_kegg(model, manifest: dict) -> None:
    kegg2model, src = km.load_reaction_crosswalk(model)
    c2bigg = km.load_metabolite_crosswalk()
    cpd_names = km.load_compound_names()
    rxn_names = km.build_reaction_names(model)
    print(f"crosswalk: {src['total']} KEGG reactions -> iCHO3K "
          f"(+{src['mnx_added']} via MetaNetX); {len(cpd_names)} compound names")

    pathways = []
    for line in _download(KEGG_LIST_URL).decode().splitlines():
        pid, _, title = line.partition("\t")
        pid = pid.replace("path:", "")
        m = re.match(r"[a-z]+(\d+)$", pid)
        if not m or int(m.group(1)) >= 1000:   # skip global/overview maps (01xxx)
            continue
        title = re.sub(r"\s*-\s*Homo sapiens.*$", "", title)
        pathways.append((pid, int(m.group(1)), title))

    kept = 0
    for pid, num, title in pathways:
        try:
            kgml = km.fetch_kgml(pid)
            emap, stats = km.build_map_from_kgml(pid, kgml, kegg2model, c2bigg,
                                                 cpd_names, rxn_names)
        except Exception as exc:
            print(f"  SKIP {pid}: {str(exc)[:50]}")
            continue
        if stats.n_reactions_mapped < MIN_MAPPED:
            continue
        out = MAPS_DIR / f"kegg_{pid}.json"
        out.write_text(json.dumps(emap))
        tier = "Medium" if stats.coverage >= MEDIUM_COVERAGE else "Low"
        amb = stats.n_ambiguous / max(stats.n_reactions_mapped, 1)
        manifest.setdefault(kegg_category(num), []).append({
            "pathway": stats.title, "file": out.name, "source": f"KEGG {pid} (human) × iCHO3K",
            "tier": tier, "coverage": round(stats.coverage, 3),
            "n_reactions": stats.n_reactions_total, "n_mapped": stats.n_reactions_mapped,
            "ambiguous_frac": round(amb, 3),
        })
        kept += 1
        time.sleep(0.15)   # be polite to KEGG
    print(f"KEGG maps kept: {kept} (of {len(pathways)} human pathways attempted)")


def main() -> None:
    model = model_setup.get_base_model()
    model_rxn = {r.id for r in model.reactions}
    MAPS_DIR.mkdir(parents=True, exist_ok=True)
    vendor_escher_js()

    manifest: dict[str, list] = {}
    fetch_recon1(model_rxn, manifest)
    fetch_kegg(model, manifest)

    # order maps within each category by reliability then coverage
    for items in manifest.values():
        items.sort(key=lambda m: (m["tier"] != "High", m["tier"] != "Medium", -m["coverage"]))

    (MAPS_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2))
    total = sum(len(v) for v in manifest.values())
    tiers = {"High": 0, "Medium": 0, "Low": 0}
    for items in manifest.values():
        for it in items:
            tiers[it["tier"]] += 1
    print(f"\nWrote {total} maps across {len(manifest)} categories -> {MAPS_DIR}")
    print(f"Reliability: {tiers['High']} High · {tiers['Medium']} Medium · {tiers['Low']} Low")


if __name__ == "__main__":
    main()
