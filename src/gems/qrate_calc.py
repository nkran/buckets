"""Compute specific exchange (q-)rates from bioprocess data and turn them into
per-sample exchange-reaction bounds for the GSMM.

Specific rate over an interval [t0, t1]:

    q_S = (S1 - S0) / IVCD           with   IVCD = ∫ VCD dt ≈ trapezoid

IVCD (integral of viable cell density) is the biomass-time the metabolite change
is normalised against. With cobra's sign convention:
  * a metabolite that is *consumed*  (S decreasing) -> uptake  -> NEGATIVE flux
  * a metabolite that is *produced*  (S increasing) -> secretion -> POSITIVE flux

Because ``q_S`` already carries the right sign (ΔS/IVCD), no extra negation is
needed: consumed metabolite gives ΔS<0 → negative q → uptake.

Units: VCD in 1e6 cells/mL, metabolites in mM. The absolute scale is arbitrary for
a generic GSMM (biomass is dimensionless growth), so rates are additionally rescaled
to a model-friendly magnitude via ``FLUX_SCALE``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C

# Rescale raw mM / (1e6 cell·mL·day) rates into a magnitude the GSMM tolerates
# (order ~mmol/gDW/h). Placeholder calibration so measured exchanges are non-trivial
# against the defined medium; does not change the relative uptake/secretion pattern.
# TODO(real-data): calibrate from an actual cell dry-weight / volume conversion so
# rates are true mmol/gDW/h, and tighten the medium so exchanges become binding.
FLUX_SCALE = 30.0   # TODO(real-data): placeholder; ~0.119 for 350 pg DW/cell (README section 8, step 5)
# fractional tolerance band applied around each point estimate to form lb/ub
QRATE_TOL = 0.20


def compute_qrates(metadata: pd.DataFrame) -> pd.DataFrame:
    """Return long-format q-rates per sample interval.

    Each vessel = (project, clone, condition). Rates are computed between
    consecutive timepoints and attributed to the *later* timepoint sample
    (the interval ending at that sample), consistent with a soft-sensor that
    reads the state accumulated up to a day.

    Columns: project, clone, condition, timepoint, metabolite, q, qP, ivcd.
    """
    out = []
    group_keys = ["project", "clone", "condition"]
    for (proj, clone, cond), g in metadata.groupby(group_keys):
        g = g.sort_values("timepoint")
        days = g["timepoint"].to_numpy(dtype=float)
        vcd = g["VCD"].to_numpy(dtype=float)
        for i in range(1, len(g)):
            dt = days[i] - days[i - 1]
            if dt <= 0:
                continue
            ivcd = 0.5 * (vcd[i] + vcd[i - 1]) * dt  # trapezoid biomass-time
            if ivcd <= 0:
                continue
            tp = g.iloc[i]["timepoint"]
            rec = {"project": proj, "clone": clone, "condition": cond,
                   "timepoint": tp, "ivcd": ivcd}
            # metabolite specific rates
            for met in C.METABOLITE_EXCHANGE:
                s0 = g.iloc[i - 1][met]
                s1 = g.iloc[i][met]
                q = (s1 - s0) / ivcd * FLUX_SCALE
                out.append({**rec, "metabolite": met, "q": q})
            # specific productivity qP from titer (mg/L per biomass-time)
            p0 = g.iloc[i - 1]["Titer"]
            p1 = g.iloc[i]["Titer"]
            qP = (p1 - p0) / ivcd * FLUX_SCALE
            for r in out[-len(C.METABOLITE_EXCHANGE):]:
                r["qP"] = qP
    df = pd.DataFrame(out)
    return df


def qrate_bounds_for_sample(qrates: pd.DataFrame, project, clone, condition,
                            timepoint, tol: float = QRATE_TOL) -> dict[str, tuple[float, float]]:
    """Return {exchange_reaction_id: (lb, ub)} for one sample.

    A tolerance band around the point estimate keeps the LP feasible while still
    constraining direction and magnitude. Uptake (q<0): bounds are both negative;
    secretion (q>0): both positive; the band always brackets the estimate.
    """
    sel = qrates[
        (qrates["project"] == project) & (qrates["clone"] == clone)
        & (qrates["condition"] == condition) & (qrates["timepoint"] == timepoint)
    ]
    bounds: dict[str, tuple[float, float]] = {}
    for _, row in sel.iterrows():
        rid = C.METABOLITE_EXCHANGE.get(row["metabolite"])
        if rid is None:
            continue
        q = float(row["q"])
        band = abs(q) * tol
        lb, ub = q - band, q + band
        if lb > ub:
            lb, ub = ub, lb
        bounds[rid] = (lb, ub)
    return bounds


# ---------------------------------------------------------------------------
# Measured rates supplied directly (already mmol/gDW/h, keyed by exchange reaction)
# ---------------------------------------------------------------------------
UNMEASURED_UPTAKE_CAP = 0.05   # mmol/gDW/h allowed for medium components with no measurement
RECONCILE_BAND = 0.02          # relative band kept around each reconciled rate


def reconcile_rates(model, measured: dict[str, float], *, min_growth: float = 0.0,
                    biomass_rxn: str = C.BIOMASS_RXN,
                    unmeasured_cap: float = UNMEASURED_UPTAKE_CAP,
                    band: float = RECONCILE_BAND):
    """Turn measured exchange rates into bounds the model can actually satisfy.

    Raw measurements are noisy (e.g. an essential amino acid measured as slightly
    secreted), so fixing every exchange at its measured value is often infeasible.
    This solves one LP for the exchange profile **closest to the measurements**
    (minimum total absolute deviation) that still obeys steady state and supports
    growth of at least ``min_growth``. Medium components without a measurement are
    capped at a small trace uptake so they cannot act as hidden carbon sources.

    Returns ``(bounds, report)``: ``bounds`` maps exchange id -> (lb, ub) for
    :func:`gems.gimme_fba.run_gimme`; ``report`` is a DataFrame with the measured
    value, the value used, and the deviation for every measured exchange. Returns
    ``(None, None)`` if no feasible profile exists.
    """
    from optlang.symbolics import Zero

    rxn_ids = {r.id for r in model.reactions}
    measured = {k: float(v) for k, v in measured.items() if k in rxn_ids and v == v}
    free = set(C.MEDIUM_FREE)
    caps: dict[str, tuple[float, float]] = {}

    with model:
        for ex in model.exchanges:
            if ex.id not in measured and ex.id not in free and ex.lower_bound < -unmeasured_cap:
                caps[ex.id] = (-unmeasured_cap, ex.upper_bound)
                ex.lower_bound = -unmeasured_cap
        model.reactions.get_by_id(biomass_rxn).lower_bound = max(min_growth, 0.0)

        slack = []
        for rid, q in measured.items():
            rxn = model.reactions.get_by_id(rid)
            rxn.bounds = (-1000.0, 1000.0)
            pos = model.problem.Variable(f"dev_pos_{rid}", lb=0)
            neg = model.problem.Variable(f"dev_neg_{rid}", lb=0)
            cons = model.problem.Constraint(rxn.flux_expression - pos + neg, lb=q, ub=q,
                                            name=f"fit_{rid}")
            model.add_cons_vars([pos, neg, cons])
            slack += [pos, neg]
        model.objective = model.problem.Objective(Zero, direction="min")
        model.objective.set_linear_coefficients({v: 1.0 for v in slack})
        sol = model.optimize()
        if sol.status != "optimal":
            return None, None
        used = {rid: float(sol.fluxes[rid]) for rid in measured}

    bounds = dict(caps)
    rows = []
    for rid, q in measured.items():
        u = used[rid]
        w = max(band * abs(u), 1e-4)
        bounds[rid] = (u - w, u + w)
        rows.append({"reaction_id": rid, "measured": q, "used": u, "deviation": u - q})
    return bounds, pd.DataFrame(rows)
