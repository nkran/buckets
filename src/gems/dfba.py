"""Dynamic FBA of a fed-batch bioreactor, for comparing FB and intensified FB (iFB).

Every time step the genome-scale model sees the current medium concentrations and
answers: given what can be taken up right now, how fast can the cells grow, and what do
they take up and secrete? Concentrations and cell density are then advanced and the
next step is solved. Feeds are added as daily boluses.

The per-cell metabolism is the same in both processes; FB and iFB differ only in
seeding density and feeding. Differences in the outcome therefore come from the
environment the cells create: at high density nutrients run out sooner and by-products
build up faster.

Assumptions (all adjustable in :class:`Physiology` or the medium table):
  * uptake rate of each medium component = max_uptake * S / (Km + S), and never more
    than what is left in the reactor during the step;
  * glucose and glutamine are taken up at that rate whether or not growth needs it
    (CHO cells do not throttle them), so surplus carbon and nitrogen leave as lactate,
    ammonia and alanine, as in real cultures;
  * only typical CHO by-products may be secreted (``SECRETABLE``);
  * lactate and ammonia slow metabolism (non-competitive inhibition of all uptakes)
    and raise the death rate; cells that cannot grow die faster;
  * uptake capacity declines with culture age (``capacity_halflife_days``), the usual
    empirical stand-in for the metabolic slow-down of stationary-phase CHO cells;
  * vitamins and other trace components are never limiting.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C
from . import model_setup
from .gimme_fba import cold_start

MEDIA_DIR = C.ROOT / "media"
DEFAULT_MEDIUM = MEDIA_DIR / "generic_cd_medium.csv"
MEDIUM_COLUMNS = ["exchange_id", "name", "basal_mM", "feed_mM", "max_uptake", "km_mM",
                  "forced_uptake"]

# by-products tracked in the reactor even if not in the medium
BYPRODUCTS = {"EX_lac_L_e": "lactate", "EX_nh4_e": "ammonium", "EX_ala_L_e": "alanine",
              "EX_glu_L_e": "glutamate", "EX_gly_e": "glycine", "EX_pyr_e": "pyruvate"}
# typical CHO secretions; every other exchange may only take up (or is closed)
SECRETABLE = ["EX_lac_L_e", "EX_co2_e", "EX_hco3_e", "EX_nh4_e", "EX_ala_L_e", "EX_glu_L_e",
              "EX_gly_e", "EX_pyr_e", "EX_for_e", "EX_h2o_e", "EX_h_e", "EX_asp_L_e",
              "EX_pro_L_e", "EX_5oxpro_e", "EX_cit_e"]
TRACE_UPTAKE = 0.002       # mmol/gDW/h allowed for vitamins / cofactors (non-limiting)
PARSIMONY = 1e-5           # weight of total flux in the objective (growth still dominates)


@dataclass
class Process:
    name: str
    seed_vcd: float               # 1e6 viable cells / mL at day 0
    feed_start_day: float         # first daily bolus feed
    feed_pct_per_day: float       # bolus volume, % of current volume per day
    glucose_setpoint_mM: float    # daily top-up to this level (0 = off)
    duration_days: float = 14.0


@dataclass
class Physiology:
    cell_dw_pg: float = 350.0          # dry weight per cell
    uptake_scale: float = 1.0          # multiplies every max_uptake in the medium table
    o2_max: float = 1.5                # mmol/gDW/h, respiratory capacity
    atp_maintenance: float = 1.0       # mmol ATP/gDW/h spent regardless of growth
    death_rate: float = 0.002          # 1/h baseline
    starvation_death: float = 0.01     # extra 1/h while the cells cannot grow
    ki_lactate_mM: float = 60.0        # lactate level halving metabolism
    ki_ammonia_mM: float = 15.0        # ammonium level halving metabolism
    lactate_uptake_max: float = 0.1    # mmol/gDW/h, lactate re-consumption capacity
    lactate_km_mM: float = 5.0
    capacity_halflife_days: float = 4.0  # uptake capacity halves every N days of culture (0 = off)
    dt_h: float = 4.0                  # time step (must divide 24)


def default_processes() -> dict[str, Process]:
    return {
        "FB": Process("FB", seed_vcd=0.5, feed_start_day=3, feed_pct_per_day=4,
                      glucose_setpoint_mM=15),
        "iFB": Process("iFB", seed_vcd=5.0, feed_start_day=1, feed_pct_per_day=6,
                       glucose_setpoint_mM=15),
    }


def load_medium(path: str | Path = DEFAULT_MEDIUM) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = [c for c in MEDIUM_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"medium table is missing columns {missing}")
    df["forced_uptake"] = df["forced_uptake"].astype(str).str.lower().isin(["true", "1", "yes"])
    return df[MEDIUM_COLUMNS]


# ---------------------------------------------------------------------------
# simulation
# ---------------------------------------------------------------------------
def _close_medium(model, phys: Physiology, medium_ids: set[str]) -> None:
    """Exchanges: closed unless medium, trace nutrient, free, or allowed secretion."""
    for r in model.exchanges:
        r.bounds = (0.0, 0.0)
    for rid in C.BASE_MEDIUM_UPTAKES:
        if rid in model.reactions and rid not in medium_ids:
            model.reactions.get_by_id(rid).lower_bound = -TRACE_UPTAKE
    for rid in C.MEDIUM_FREE:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds = (-1000.0, 1000.0)
    for rid in SECRETABLE:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).upper_bound = 1000.0
    model.reactions.get_by_id("EX_o2_e").bounds = (-phys.o2_max, 0.0)
    if "ATPM" in model.reactions:
        model.reactions.get_by_id("ATPM").lower_bound = phys.atp_maintenance


def _set_objective(model, biomass_rxn: str = C.BIOMASS_RXN) -> None:
    """Maximise growth minus a tiny cost on total flux, in one LP.

    Plain FBA leaves every non-limiting uptake arbitrary, so the solver may take far
    more amino acid than growth needs and burn it to ammonia. The small flux cost makes
    the cells take only what they need (forced glucose/glutamine excepted) and picks
    the most economical by-product route, without trading off any meaningful growth.
    """
    from optlang.symbolics import Zero
    model.objective = model.problem.Objective(Zero, direction="max")
    coefs = {}
    for r in model.reactions:
        coefs[r.forward_variable] = -PARSIMONY
        coefs[r.reverse_variable] = -PARSIMONY
    bio = model.reactions.get_by_id(biomass_rxn)
    coefs[bio.forward_variable] = 1.0 - PARSIMONY
    coefs[bio.reverse_variable] = -1.0 - PARSIMONY
    model.objective.set_linear_coefficients(coefs)


def simulate(medium: pd.DataFrame, process: Process, phys: Physiology | None = None,
             model=None) -> dict:
    """Run one bioreactor simulation. Returns {"timecourse": DataFrame, "bounds": [...]}.

    ``bounds[i]`` holds the exchange bounds used at row ``i`` of the timecourse, so any
    time point can be re-solved later for full internal fluxes (:func:`fluxes_at`).
    """
    phys = phys or Physiology()
    if 24 % phys.dt_h:
        raise ValueError("dt_h must divide 24")
    model = model if model is not None else model_setup.get_base_model()
    medium = medium.copy()
    medium = medium[medium["exchange_id"].isin({r.id for r in model.reactions})]
    med = {r.exchange_id: r for r in medium.itertuples(index=False)}

    species = list(dict.fromkeys(list(med) + list(BYPRODUCTS)))
    S = {s: float(med[s].basal_mM) if s in med else 0.0 for s in species}
    dw = phys.cell_dw_pg * 1e-3                       # gDW/L per 1e6 cells/mL
    X, Xd, ivcd = process.seed_vcd * dw, 0.0, 0.0
    dt = phys.dt_h
    n_steps = int(round(process.duration_days * 24 / dt))
    biomass = model.reactions.get_by_id(C.BIOMASS_RXN)
    rows, bounds_log = [], []

    cold_start(model)
    with model:
        _close_medium(model, phys, set(med))
        _set_objective(model)
        for step in range(n_steps + 1):
            t_h = step * dt
            fed = 0.0
            # ---- daily bolus feed + glucose top-up (at day boundaries) -------------
            if step > 0 and t_h % 24 == 0 and t_h >= process.feed_start_day * 24:
                f = process.feed_pct_per_day / 100.0
                if f > 0:
                    for s in species:
                        feed_c = float(med[s].feed_mM) if s in med else 0.0
                        S[s] = (S[s] + f * feed_c) / (1 + f)
                    X, Xd = X / (1 + f), Xd / (1 + f)
                    fed = f
                if process.glucose_setpoint_mM > 0 and "EX_glc_e" in S \
                        and S["EX_glc_e"] < process.glucose_setpoint_mM:
                    S["EX_glc_e"] = process.glucose_setpoint_mM

            # ---- uptake limits for this step ---------------------------------------
            inhib = (phys.ki_lactate_mM / (phys.ki_lactate_mM + S["EX_lac_L_e"])
                     * phys.ki_ammonia_mM / (phys.ki_ammonia_mM + S["EX_nh4_e"]))
            if phys.capacity_halflife_days > 0:
                inhib *= 0.5 ** ((t_h / 24) / phys.capacity_halflife_days)
            step_bounds: dict[str, tuple[float, float]] = {}
            for rid, row in med.items():
                s = S[rid]
                cap = row.max_uptake * phys.uptake_scale * inhib * s / (row.km_mM + s)
                if X > 0:
                    cap = min(cap, s / (X * dt))
                ub = 1000.0 if rid in SECRETABLE else 0.0
                step_bounds[rid] = (-cap, -cap) if row.forced_uptake else (-cap, ub)
            if "EX_lac_L_e" not in med:
                lac = S["EX_lac_L_e"]
                cap = phys.lactate_uptake_max * inhib * lac / (phys.lactate_km_mM + lac)
                if X > 0:
                    cap = min(cap, lac / (X * dt))
                step_bounds["EX_lac_L_e"] = (-cap, 1000.0)
            for rid, b in step_bounds.items():
                model.reactions.get_by_id(rid).bounds = b

            model.slim_optimize()
            relaxed = False
            if model.solver.status != "optimal":   # cannot use all forced uptake: relax it
                for rid, row in med.items():
                    if row.forced_uptake:
                        lb = step_bounds[rid][0]
                        step_bounds[rid] = (lb, 0.0)
                        model.reactions.get_by_id(rid).bounds = (lb, 0.0)
                relaxed = True
                model.slim_optimize()
            if model.solver.status == "optimal":
                mu = max(biomass.flux, 0.0)
                v = {s: model.reactions.get_by_id(s).flux for s in species}
                # growth-limiting nutrient: at its uptake cap with the largest shadow price
                limiting, best = "", 0.0
                for rid, row in med.items():
                    lb = step_bounds[rid][0]
                    if lb < 0 and not row.forced_uptake and abs(v[rid] - lb) < 1e-9:
                        rc = abs(model.reactions.get_by_id(rid).reduced_cost)
                        if rc > best + 1e-6:
                            limiting, best = row.name, rc
                if not limiting and mu < 1e-4:
                    limiting = "(not growing)"
            else:
                mu, v, limiting = 0.0, {s: 0.0 for s in species}, "(no feasible metabolism)"

            vcd = X / dw
            rows.append({"day": t_h / 24, "vcd": vcd, "viability": X / (X + Xd) if X + Xd else 0,
                         "growth_rate": mu, "ivcd": ivcd, "limiting": limiting,
                         "inhibition": inhib, "feed_added": fed, "forced_relaxed": relaxed,
                         **{f"conc:{s}": S[s] for s in species},
                         **{f"rate:{s}": v[s] for s in species}})
            bounds_log.append(step_bounds)
            if step == n_steps:
                break

            # ---- integrate over the step -------------------------------------------
            kd = phys.death_rate * (1 + S["EX_lac_L_e"] / phys.ki_lactate_mM
                                    + S["EX_nh4_e"] / phys.ki_ammonia_mM)
            if mu < 1e-3:
                kd += phys.starvation_death
            g = mu - kd
            integral = X * dt if abs(g) < 1e-12 else X * (math.exp(g * dt) - 1) / g
            for s in species:
                S[s] = max(S[s] + v[s] * integral, 0.0)
            Xd += kd * integral
            ivcd += (integral / dw) / 24.0            # 1e6 cells·day/mL
            X = X * math.exp(g * dt)

    tc = pd.DataFrame(rows)
    return {"timecourse": tc, "bounds": bounds_log, "process": asdict(process),
            "physiology": asdict(phys)}


def fluxes_at(result: dict, day: float, model=None) -> pd.Series:
    """Full internal flux distribution at the time point nearest ``day`` (pFBA:
    maximum growth, then the most economical flux pattern)."""
    from cobra.flux_analysis import pfba

    model = model if model is not None else model_setup.get_base_model()
    tc = result["timecourse"]
    i = int((tc["day"] - day).abs().idxmin())
    phys = Physiology(**result["physiology"])
    cold_start(model)
    with model:
        _close_medium(model, phys, set(result["bounds"][i]))
        for rid, b in result["bounds"][i].items():
            model.reactions.get_by_id(rid).bounds = b
        model.objective = C.BIOMASS_RXN
        try:
            return pfba(model).fluxes
        except Exception:
            return pd.Series(dtype=float)


def summarise(result: dict) -> dict:
    """Headline numbers for one run."""
    tc = result["timecourse"]
    lim = tc.loc[tc["limiting"] != "", ["day", "limiting"]]
    first = lim.iloc[0] if len(lim) else None
    return {
        "peak VCD (1e6/mL)": tc["vcd"].max(),
        "day of peak": tc.loc[tc["vcd"].idxmax(), "day"],
        "final IVCD (1e6 cell·d/mL)": tc["ivcd"].iloc[-1],
        "final viability": tc["viability"].iloc[-1],
        "peak lactate (mM)": tc["conc:EX_lac_L_e"].max(),
        "final ammonium (mM)": tc["conc:EX_nh4_e"].iloc[-1],
        "first limitation": f"{first['limiting']} (day {first['day']:.1f})" if first is not None else "none",
    }


def depletion_days(result: dict, medium: pd.DataFrame, threshold_frac: float = 0.05) -> pd.Series:
    """First day each medium component falls below ``threshold_frac`` of its basal level."""
    tc = result["timecourse"]
    out = {}
    for r in medium.itertuples(index=False):
        col = f"conc:{r.exchange_id}"
        if col not in tc or r.basal_mM <= 0:
            continue
        low = tc.loc[tc[col] < threshold_frac * r.basal_mM, "day"]
        out[r.name] = float(low.iloc[0]) if len(low) else np.nan
    return pd.Series(out, name="depleted on day").sort_values()
