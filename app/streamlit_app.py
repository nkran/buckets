"""Interactive Escher explorer for iCHO3K flux + expression.

Run:  streamlit run app/streamlit_app.py

Browse pathway maps and overlay, per clone/timepoint, either predicted FLUX or gene
EXPRESSION on the reaction arrows (color + thickness), or compare them side by side.
Reads the matrices exported by the pipeline (results/fluxes/) — no cobra/solver needed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from gems import config as C
from gems import escher_view as ev

MAPS_DIR = C.ROOT / "escher_maps"
META_CANDIDATES = [C.DATA / "raw" / "metadata.csv", C.DATA / "synthetic" / "metadata.csv"]

st.set_page_config(page_title="iCHO3K Flux/Expression Explorer", layout="wide")


# ---------------------------------------------------------------------------
# Loaders (cached)
# ---------------------------------------------------------------------------
def find_datasets() -> dict[str, Path]:
    """Every pipeline run under results/ that has explorer matrices: name -> fluxes dir."""
    found = {}
    for p in sorted(C.RESULTS.glob("**/flux_matrix.parquet")):
        d = p.parent
        name = "default (results/)" if d.parent == C.RESULTS else str(d.parent.relative_to(C.RESULTS))
        found[name] = d
    return found


@st.cache_data(show_spinner=False)
def load_matrix(kind: str, flux_dir: str) -> pd.DataFrame | None:
    for ext in ("parquet", "csv"):
        p = Path(flux_dir) / f"{kind}.{ext}"
        if p.exists():
            df = pd.read_parquet(p) if ext == "parquet" else pd.read_csv(p, index_col=0)
            df.index = df.index.astype(str)
            return df
    return None


@st.cache_data(show_spinner=False)
def load_manifest() -> dict:
    p = MAPS_DIR / "manifest.json"
    return json.loads(p.read_text()) if p.exists() else {}


@st.cache_data(show_spinner=False)
def load_metadata(sample_ids: tuple[str, ...], flux_dir: str) -> pd.DataFrame | None:
    want = set(sample_ids)
    summary = Path(flux_dir) / "sample_summary.csv"      # written by the pipeline run
    if summary.exists():
        m = pd.read_csv(summary, index_col=0)
        m.index = m.index.astype(str)
        if want & set(m.index):
            return m
    for p in META_CANDIDATES:
        if p.exists():
            m = pd.read_csv(p, index_col=0)
            m.index = m.index.astype(str)
            if want & set(m.index):
                return m
    return None


@st.cache_data(show_spinner=True)
def render_html(map_file: str, sample: str, mode: str, kind: str, flux_dir: str) -> tuple[str, float]:
    matrix = load_matrix(kind, flux_dir)
    emap = json.loads((MAPS_DIR / map_file).read_text())
    data = ev.reaction_data_for_sample(matrix, sample)
    cov = ev.coverage_for_map(emap, ev.resolve_aliases(MAPS_DIR / map_file, data))
    html = ev.build_escher_html(MAPS_DIR / map_file, data, mode=mode)
    return html, cov


def sample_label(sid: str, meta: pd.DataFrame | None) -> str:
    if meta is None or sid not in meta.index:
        return sid
    r = meta.loc[sid]
    bits = [str(r.get(k)) for k in ("project", "clone", "condition", "timepoint") if k in r]
    return f"{sid}  ({' · '.join(bits)})"


# ---------------------------------------------------------------------------
# App
# ---------------------------------------------------------------------------
st.title("🧬 iCHO3K Flux / Expression Explorer")

datasets = find_datasets()
manifest = load_manifest()
if datasets:
    with st.sidebar:
        st.header("Controls")
        ds_name = st.selectbox("Dataset", list(datasets), help="one entry per pipeline run under results/")
    FLUX_DIR = str(datasets[ds_name])
    flux = load_matrix("flux_matrix", FLUX_DIR)
    expr = load_matrix("rxn_expression", FLUX_DIR)
else:
    flux = expr = None

if flux is None or not manifest:
    st.error(
        "Missing inputs. Run the pipeline (to create `results/fluxes/flux_matrix.parquet` "
        "and `rxn_expression.parquet`) and `python scripts/prepare_escher_maps.py` "
        "(to create `escher_maps/`), then reload."
    )
    st.stop()

samples = list(flux.columns)
meta = load_metadata(tuple(samples), FLUX_DIR)

with st.sidebar:
    st.caption(f"{flux.shape[0]} reactions · {len(samples)} samples")

    sample = st.selectbox("Sample", samples, format_func=lambda s: sample_label(s, meta))

    # reliability filter
    tiers = st.multiselect("Reliability", ["High", "Medium", "Low"],
                           default=["High", "Medium", "Low"],
                           help="High = curated RECON1 map · Medium/Low = KEGG-derived "
                                "(by fraction of pathway reactions mapped to iCHO3K)")

    # pathway map grouped by category, filtered by tier
    _badge = {"High": "🟢", "Medium": "🟡", "Low": "🔴"}
    labels, files, entries = [], [], []
    for category in sorted(manifest):
        for it in manifest[category]:
            if it.get("tier", "Low") not in tiers:
                continue
            labels.append(f"{_badge.get(it.get('tier'),'')} {category} — {it['pathway']} "
                          f"[{it['coverage']:.0%}]")
            files.append(it["file"])
            entries.append({**it, "category": category})
    if not files:
        st.warning("No maps match the selected reliability filter.")
        st.stop()
    idx = st.selectbox("Pathway map", range(len(labels)), format_func=lambda i: labels[i])
    map_file = files[idx]
    entry = entries[idx]

    metric = st.radio("Overlay", ["Flux", "Expression", "Compare"], index=0)
    st.markdown(
        "- **Flux**: thickness = |flux|, colour = direction "
        "(<span style='color:#2166ac'>blue −</span> / "
        "<span style='color:#b2182b'>red +</span>)\n"
        "- **Expression**: thickness + colour = reaction expression score\n"
        "- Grey = no data for that reaction on this map",
        unsafe_allow_html=True,
    )
    if expr is None and metric != "Flux":
        st.warning("No expression matrix found — only Flux is available.")
        metric = "Flux"

st.subheader(f"{entry['category']} — {entry['pathway']}")
st.caption(f"Sample: **{sample_label(sample, meta)}**")

# reliability panel
rc1, rc2, rc3, rc4 = st.columns(4)
rc1.metric("Reliability", f"{_badge.get(entry.get('tier'),'')} {entry.get('tier','?')}")
rc2.metric("Pathway coverage", f"{entry['coverage']:.0%}",
           help=f"{entry['n_mapped']} of {entry['n_reactions']} pathway reactions map to iCHO3K")
rc3.metric("Reactions mapped", entry["n_mapped"])
rc4.metric("Source", "RECON1" if "RECON1" in entry["source"] else "KEGG")
with st.expander("How to read this map's reliability"):
    st.markdown(
        f"- **Source**: {entry['source']}\n"
        f"- **Pathway coverage**: {entry['n_mapped']}/{entry['n_reactions']} "
        f"({entry['coverage']:.0%}) of the pathway's reactions exist in iCHO3K and carry an "
        f"overlay; the rest render grey.\n"
        f"- **Ambiguous mappings**: {entry.get('ambiguous_frac', 0):.0%} of mapped reactions "
        f"had several iCHO3K matches (one chosen) — mainly compartment variants.\n"
        "- **High** tier = curated RECON1 human map (direct BiGG ids, hand-drawn). "
        "**Medium/Low** = KEGG layout bridged to iCHO3K via BiGG; trust the connected "
        "*core* of the pathway, treat sparse peripheral reactions cautiously."
    )

if metric == "Compare":
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Flux**")
        html, cov = render_html(map_file, sample, "flux", "flux_matrix", FLUX_DIR)
        st.caption(f"map coverage: {cov:.0%}")
        components.html(html, height=560, scrolling=False)
    with c2:
        st.markdown("**Expression**")
        html, cov = render_html(map_file, sample, "expression", "rxn_expression", FLUX_DIR)
        st.caption(f"map coverage: {cov:.0%}")
        components.html(html, height=560, scrolling=False)
else:
    kind = "flux_matrix" if metric == "Flux" else "rxn_expression"
    mode = "flux" if metric == "Flux" else "expression"
    html, cov = render_html(map_file, sample, mode, kind, FLUX_DIR)
    st.caption(f"map coverage: {cov:.0%}  ·  drag to pan, scroll to zoom")
    components.html(html, height=680, scrolling=False)

    # top reactions on this map for the selected sample
    matrix = flux if metric == "Flux" else expr
    emap = json.loads((MAPS_DIR / map_file).read_text())
    rids = [r["bigg_id"] for r in emap[1]["reactions"].values()]
    col = matrix[sample].reindex(rids).dropna()
    if not col.empty:
        top = col.reindex(col.abs().sort_values(ascending=False).index).head(15)
        with st.expander(f"Top reactions by |{metric.lower()}| on this map"):
            st.dataframe(top.rename(metric).to_frame(), width="stretch")
