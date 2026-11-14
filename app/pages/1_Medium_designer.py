"""Medium & process designer: change the medium and compare FB vs iFB (dynamic FBA).

Opened as the second page of:  streamlit run app/streamlit_app.py
"""

from __future__ import annotations

import io
import json
import threading
from dataclasses import asdict, fields

import altair as alt
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from gems import config as C
from gems import dfba
from gems import escher_view as ev
from gems import model_setup

st.set_page_config(page_title="Medium designer", layout="wide")
COLORS = {"FB": "#4a6fa5", "iFB": "#d1782f"}
MAPS_DIR = C.ROOT / "escher_maps"


# ---------------------------------------------------------------------------
# cached helpers
# ---------------------------------------------------------------------------
@st.cache_resource(show_spinner=False)
def get_model():
    return model_setup.get_base_model()


@st.cache_resource(show_spinner=False)
def model_lock():
    # one shared model object: serialise solves if two browser tabs run at once
    return threading.Lock()


@st.cache_resource(show_spinner=False)
def exchange_catalogue() -> pd.DataFrame:
    m = get_model()
    rows = [{"exchange_id": r.id, "name": (r.reactants[0].name if r.reactants else r.id)}
            for r in m.exchanges]
    return pd.DataFrame(rows).sort_values("name")


@st.cache_data(show_spinner=False, max_entries=40)
def run_one(medium_csv: str, proc_json: str, phys_json: str) -> dict:
    """Simulate one process (cached per process, so editing iFB does not re-run FB)."""
    medium = dfba.load_medium(io.StringIO(medium_csv))
    with model_lock():
        return dfba.simulate(medium, dfba.Process(**json.loads(proc_json)),
                             dfba.Physiology(**json.loads(phys_json)), model=get_model())


@st.cache_data(show_spinner=False, max_entries=60)
def snapshot(run_key: str, proc: str, day: float, _result: dict) -> pd.Series:
    with model_lock():
        return dfba.fluxes_at(_result, day, model=get_model())


def medium_presets() -> dict[str, str]:
    return {p.stem.replace("_", " "): str(p) for p in sorted(dfba.MEDIA_DIR.glob("*.csv"))}


# ---------------------------------------------------------------------------
# page
# ---------------------------------------------------------------------------
st.title("🧪 Medium & process designer — FB vs iFB")
with st.expander("What this page does (read first)", expanded=False):
    st.markdown(
        """
Each process is simulated as a **dynamic FBA** of the bioreactor. Every few hours
the iCHO3K model sees the current medium concentrations, takes up nutrients (rate =
max uptake × S/(Km+S), never more than is left), grows as fast as those nutrients allow,
and secretes by-products. Concentrations and cell density are then advanced and feeds
added as daily boluses.

**FB vs iFB.** The cells' metabolism is identical in both; the processes differ only in
seeding density and feeding (sidebar). Any difference in the outcome comes from the
environment the cells create: at high density nutrients run out sooner and lactate and
ammonia build up faster.

**How to use it.** Edit the medium table, press **Run**, then **Pin as baseline** and
change something — every chart then shows the new run against the pinned one.

**How far to trust it.** Default uptake limits come from measured CHO fed-batch rates
(the ZeLa dataset); concentrations of the default medium are illustrative and growth
inhibition/death parameters are rough. Absolute numbers need calibrating against your
own runs (*Advanced physiology* in the sidebar); the *direction* of the response to a
medium change is the useful output. Product (antibody) formation is not modelled —
IVCD is the proxy.
"""
    )

# ---- sidebar: medium source, processes, physiology ------------------------
with st.sidebar:
    st.header("Medium")
    presets = medium_presets()
    preset = st.selectbox("Start from", list(presets),
                          help="every CSV in the media/ folder appears here")
    upload = st.file_uploader("…or upload a medium CSV", type="csv")
    src_key = f"upload:{upload.name}" if upload else f"preset:{preset}"
    if st.session_state.get("medium_src") != src_key:
        st.session_state["medium"] = dfba.load_medium(upload if upload else presets[preset])
        st.session_state["medium_src"] = src_key
        st.session_state["editor_v"] = st.session_state.get("editor_v", 0) + 1

    st.header("Processes")
    defaults = dfba.default_processes()
    duration = st.slider("Culture length (days)", 6, 21, 14)
    procs = {}
    for name, p in defaults.items():
        with st.expander(name, expanded=False):
            procs[name] = dfba.Process(
                name=name,
                seed_vcd=st.number_input("Seeding density (10⁶ cells/mL)", 0.1, 50.0,
                                         p.seed_vcd, 0.1, key=f"{name}_seed"),
                feed_start_day=st.number_input("First feed (day)", 0.0, 10.0,
                                               float(p.feed_start_day), 1.0, key=f"{name}_fs"),
                feed_pct_per_day=st.number_input("Daily feed (% of volume)", 0.0, 20.0,
                                                 float(p.feed_pct_per_day), 0.5, key=f"{name}_fp"),
                glucose_setpoint_mM=st.number_input("Glucose top-up to (mM, 0 = off)", 0.0, 60.0,
                                                    float(p.glucose_setpoint_mM), 1.0,
                                                    key=f"{name}_gs"),
                duration_days=float(duration),
            )

    with st.expander("Advanced physiology", expanded=False):
        d = dfba.Physiology()
        help_txt = {
            "cell_dw_pg": "dry weight per cell (pg)",
            "uptake_scale": "multiplies every max uptake in the medium table",
            "o2_max": "respiratory capacity, mmol O2/gDW/h",
            "atp_maintenance": "ATP spent regardless of growth, mmol/gDW/h",
            "death_rate": "baseline death rate, 1/h",
            "starvation_death": "extra death rate while cells cannot grow, 1/h",
            "ki_lactate_mM": "lactate level that halves metabolism",
            "ki_ammonia_mM": "ammonium level that halves metabolism",
            "lactate_uptake_max": "lactate re-consumption capacity, mmol/gDW/h",
            "lactate_km_mM": "lactate uptake Km",
            "capacity_halflife_days": "uptake capacity halves every N culture days (0 = off)",
        }
        phys_vals = {}
        for f in fields(dfba.Physiology):
            if f.name == "dt_h":
                phys_vals[f.name] = float(st.select_slider("Time step (h)", [2, 3, 4, 6, 8, 12],
                                                           value=int(d.dt_h)))
            else:
                v = getattr(d, f.name)
                phys_vals[f.name] = st.number_input(f.name, value=float(v),
                                                    step=float(v) / 10 if v else 0.1,
                                                    format="%.4f", help=help_txt.get(f.name))
        phys = dfba.Physiology(**phys_vals)

# ---- 1. medium editor -----------------------------------------------------
st.subheader("1 · Medium composition")
st.caption("Basal = concentration in the starting medium; feed = concentration in the "
           "concentrated feed. Max uptake (mmol/gDW/h) and Km set how fast cells can take "
           "each component up. *Forced* = taken up at that rate even when not needed "
           "(glucose, glutamine).")
edited = st.data_editor(
    st.session_state["medium"], key=f"medium_editor_{st.session_state['editor_v']}",
    num_rows="dynamic", width="stretch", hide_index=True,
    column_config={
        "exchange_id": st.column_config.TextColumn("model exchange", disabled=True),
        "name": st.column_config.TextColumn("component", disabled=True),
        "basal_mM": st.column_config.NumberColumn("basal (mM)", min_value=0.0, format="%.2f"),
        "feed_mM": st.column_config.NumberColumn("feed (mM)", min_value=0.0, format="%.2f"),
        "max_uptake": st.column_config.NumberColumn("max uptake", min_value=0.0, format="%.4f"),
        "km_mM": st.column_config.NumberColumn("Km (mM)", min_value=0.0, format="%.2f"),
        "forced_uptake": st.column_config.CheckboxColumn("forced"),
    },
)
c1, c2, c3 = st.columns([3, 1, 1])
cat = exchange_catalogue()
cat = cat[~cat.exchange_id.isin(edited.exchange_id)]
new = c1.selectbox("Add a component", ["—"] + [f"{n}  [{i}]" for i, n in zip(cat.exchange_id, cat.name)])
if c2.button("Add", disabled=new == "—"):
    rid = new.split("[")[-1].rstrip("]")
    row = {"exchange_id": rid, "name": new.split("  [")[0], "basal_mM": 0.0, "feed_mM": 0.0,
           "max_uptake": 0.02, "km_mM": 0.1, "forced_uptake": False}
    st.session_state["medium"] = pd.concat([edited, pd.DataFrame([row])], ignore_index=True)
    st.session_state["editor_v"] += 1
    st.rerun()
c3.download_button("Download medium CSV", edited.to_csv(index=False),
                   file_name="my_medium.csv", mime="text/csv")

# ---- 2. run ---------------------------------------------------------------
st.subheader("2 · Simulate")
medium_csv = edited.dropna(subset=["exchange_id"]).to_csv(index=False)
procs_json = json.dumps({k: asdict(v) for k, v in procs.items()}, sort_keys=True)
phys_json = json.dumps(asdict(phys), sort_keys=True)
run_key = str(hash((medium_csv, procs_json, phys_json)))

b1, b2, b3 = st.columns([1, 1, 4])
if b1.button("▶ Run", type="primary"):
    out = {}
    for name, p in procs.items():
        with st.spinner(f"Simulating {name} (≈15 s)…"):
            out[name] = run_one(medium_csv, json.dumps(asdict(p), sort_keys=True), phys_json)
    st.session_state["result"] = out
    st.session_state["result_key"] = run_key
    st.session_state["result_medium"] = edited.copy()
res = st.session_state.get("result")
if res is not None and b2.button("📌 Pin as baseline"):
    st.session_state["baseline"] = res
    st.session_state["baseline_key"] = st.session_state["result_key"]
base = st.session_state.get("baseline")
if base is not None and b3.button("Clear baseline"):
    st.session_state.pop("baseline", None)
    base = None
if res is None:
    st.info("Edit the medium if you like, then press **Run**.")
    st.stop()
if st.session_state.get("result_key") != run_key:
    st.warning("Inputs changed since the last run — press **Run** to update.")
if base is not None and st.session_state.get("baseline_key") == st.session_state.get("result_key"):
    st.caption("The current run is the pinned baseline; change something and run again to compare.")
med_used = st.session_state["result_medium"]

# ---- KPIs -----------------------------------------------------------------
kpi_cols = st.columns(len(res))
for col, (name, r) in zip(kpi_cols, res.items()):
    s = dfba.summarise(r)
    sb = dfba.summarise(base[name]) if base is not None and name in base else None
    col.markdown(f"#### {name}")
    k1, k2, k3 = col.columns(3)
    for box, key, fmt in [(k1, "peak VCD (1e6/mL)", "{:.1f}"), (k2, "final IVCD (1e6 cell·d/mL)", "{:.0f}"),
                          (k3, "final viability", "{:.0%}")]:
        delta = None if sb is None else fmt.format(s[key] - sb[key])
        box.metric(key, fmt.format(s[key]), delta)
    k4, k5, k6 = col.columns(3)
    for box, key, fmt, inv in [(k4, "peak lactate (mM)", "{:.1f}", True),
                               (k5, "final ammonium (mM)", "{:.1f}", True),
                               (k6, "day of peak", "{:.1f}", False)]:
        delta = None if sb is None else fmt.format(s[key] - sb[key])
        box.metric(key, fmt.format(s[key]), delta, delta_color="inverse" if inv else "normal")
    col.caption(f"First growth limitation: **{s['first limitation']}**")


# ---- long-format helper -----------------------------------------------------
def long(results: dict, run_label: str) -> pd.DataFrame:
    frames = []
    for name, r in results.items():
        tc = r["timecourse"].copy()
        tc["viability_pct"] = 100 * tc["viability"]
        tc["process"], tc["run"] = name, run_label
        frames.append(tc)
    return pd.concat(frames, ignore_index=True)


cur = long(res, "current")
both = pd.concat([cur, long(base, "baseline")]) if base is not None else cur
names = dict(zip(med_used.exchange_id, med_used.name))
names.update({k: v for k, v in dfba.BYPRODUCTS.items() if k not in names})


def chart(col: str, title: str, unit: str):
    d = both[["day", "process", "run", col]].rename(columns={col: "value"})
    return (alt.Chart(d, title=title, height=210)
            .mark_line()
            .encode(x=alt.X("day:Q", title="day"),
                    y=alt.Y("value:Q", title=unit),
                    color=alt.Color("process:N", scale=alt.Scale(domain=list(COLORS),
                                                                 range=list(COLORS.values()))),
                    strokeDash=alt.StrokeDash("run:N", scale=alt.Scale(
                        domain=["current", "baseline"], range=[[1, 0], [5, 4]])),
                    tooltip=["process", "run", alt.Tooltip("day:Q", format=".1f"),
                             alt.Tooltip("value:Q", format=".3f")]))


tab1, tab2, tab3, tab4 = st.tabs(["Culture", "What limits growth", "Fluxes inside the cell",
                                  "Data"])

with tab1:
    main = [("vcd", "Viable cell density", "10⁶ cells/mL"), ("viability_pct", "Viability", "%"),
            ("growth_rate", "Specific growth rate", "1/h"),
            ("conc:EX_glc_e", "Glucose", "mM"), ("conc:EX_lac_L_e", "Lactate", "mM"),
            ("conc:EX_nh4_e", "Ammonium", "mM")]
    for row in (main[:3], main[3:]):
        cols = st.columns(3)
        for c, (col, t, u) in zip(cols, row):
            c.altair_chart(chart(col, t, u), use_container_width=True)
    extra = [c for c in cur.columns if c.startswith("conc:")
             and c not in {m[0] for m in main}]
    pick = st.multiselect("More components", extra,
                          default=[c for c in extra if c.endswith(("gln_L_e", "asn_L_e", "ala_L_e"))],
                          format_func=lambda c: names.get(c[5:], c[5:]))
    cols = st.columns(3)
    for i, c in enumerate(pick):
        cols[i % 3].altair_chart(chart(c, names.get(c[5:], c[5:]), "mM"), use_container_width=True)
    st.caption("Solid = current run, dashed = pinned baseline. Saw-tooth = daily bolus feeds.")

with tab2:
    st.markdown("**Which nutrient caps growth over time** (the medium component whose "
                "uptake limit binds with the largest effect on growth).")
    lim = cur[["day", "process", "limiting"]].copy()
    lim["limiting"] = lim["limiting"].replace("", "(nothing limiting)")
    st.altair_chart(
        alt.Chart(lim, height=120).mark_rect().encode(
            x=alt.X("day:Q", bin=alt.Bin(step=phys.dt_h / 24), title="day"),
            y=alt.Y("process:N", title=None),
            color=alt.Color("limiting:N", title="limiting", scale=alt.Scale(scheme="tableau20")),
            tooltip=["process", alt.Tooltip("day:Q", format=".1f"), "limiting"]),
        use_container_width=True)
    dep = pd.DataFrame({name: dfba.depletion_days(r, med_used) for name, r in res.items()})
    st.markdown("**Day each medium component falls below 5% of its starting level** "
                "(blank = never)")
    st.dataframe(dep.round(1), width="stretch")

with tab3:
    tc0 = cur[cur.process == list(res)[0]]
    day = st.slider("Day", 0.0, float(tc0.day.max()), min(5.0, float(tc0.day.max())),
                    phys.dt_h / 24)
    with st.spinner("Solving full flux distributions…"):
        snaps = {name: snapshot(st.session_state["result_key"], name, day, r)
                 for name, r in res.items()}
    model = get_model()
    sub = {r.id: (r.subsystem or "") for r in model.reactions}
    rows = []
    for name, f in snaps.items():
        if f.empty:
            continue
        for s, v in f.abs().groupby(f.index.map(sub)).sum().items():
            if s in C.PANEL_SUBSYSTEMS:
                rows.append({"process": name, "pathway": s.title(), "total |flux|": v})
    if rows:
        st.markdown(f"**Pathway activity on day {day:.1f}** (sum of |flux|, mmol/gDW/h)")
        st.altair_chart(
            alt.Chart(pd.DataFrame(rows), height=300).mark_bar().encode(
                y=alt.Y("pathway:N", sort="-x", title=None), x="total |flux|:Q",
                color=alt.Color("process:N", scale=alt.Scale(domain=list(COLORS),
                                                             range=list(COLORS.values()))),
                yOffset="process:N", tooltip=["process", "pathway",
                                              alt.Tooltip("total |flux|:Q", format=".3f")]),
            use_container_width=True)
    ex = pd.DataFrame({n: f.reindex(list(names)) for n, f in snaps.items() if not f.empty})
    ex.index = [names.get(i, i) for i in ex.index]
    st.markdown("**Uptake (−) and secretion (+) on that day**, mmol/gDW/h")
    st.dataframe(ex.round(4), width="stretch")

    manifest = MAPS_DIR / "manifest.json"
    if manifest.exists():
        items = [(f"{cat} — {it['pathway']}", it["file"])
                 for cat, its in json.loads(manifest.read_text()).items() for it in its]
        default = next((i for i, (_, f) in enumerate(items) if f == "glycolysis_tca_ppp.json"), 0)
        label = st.selectbox("Pathway map", range(len(items)), index=default,
                             format_func=lambda i: items[i][0])
        mcols = st.columns(len(snaps))
        for c, (name, f) in zip(mcols, snaps.items()):
            c.markdown(f"**{name}, day {day:.1f}** — thickness = |flux|, red/blue = direction")
            if not f.empty:
                html = ev.build_escher_html(MAPS_DIR / items[label][1],
                                            {k: float(v) for k, v in f.items()},
                                            mode="flux", height=480)
                with c:
                    components.html(html, height=490)
    else:
        st.caption("Pathway maps appear here after `python scripts/prepare_escher_maps.py`.")

with tab4:
    for name, r in res.items():
        st.download_button(f"Download {name} time course (CSV)",
                           r["timecourse"].to_csv(index=False), file_name=f"dfba_{name}.csv",
                           mime="text/csv", key=f"dl_{name}")
    st.dataframe(cur.drop(columns=[c for c in cur.columns if c.startswith("rate:")]).round(4),
                 width="stretch")
