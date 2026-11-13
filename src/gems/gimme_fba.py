"""Native COBRApy GIMME + parsimonious FBA per sample.

GIMME (Becker & Palsson, 2008), continuous variant:
  1. Score each reaction from expression (see :mod:`gene_mapping`).
  2. Choose a threshold T (a quantile of the scored reactions).
  3. Reactions scoring below T get a penalty c_i = T - score_i (data-bearing only).
  4. Required Metabolic Functionality (RMF): hold biomass ≥ fraction · biomass_max.
  5. Minimise Σ c_i·|v_i| subject to steady state + the sample's exchange bounds.

Stage 2 keeps the GIMME inconsistency at its optimum and then minimises total flux
(pFBA-style L1), yielding a unique, parsimonious, expression-consistent distribution.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from optlang.symbolics import Zero

from . import config as C


@dataclass
class GimmeResult:
    status: str
    growth: float
    inconsistency: float
    fluxes: pd.Series          # reaction_id -> flux (full model)
    growth_max: float = float("nan")   # max biomass flux the constraints allow


def _penalties(scores: dict[str, float], quantile: float) -> tuple[dict[str, float], float]:
    vals = np.array([v for v in scores.values() if not np.isnan(v)], dtype=float)
    if vals.size == 0:
        return {}, 0.0
    threshold = float(np.quantile(vals, quantile))
    pen = {rid: (threshold - s) for rid, s in scores.items()
           if not np.isnan(s) and s < threshold}
    return pen, threshold


def run_gimme(model, scores: dict[str, float],
              exchange_bounds: dict[str, tuple[float, float]] | None = None,
              quantile: float = C.GIMME_QUANTILE,
              biomass_fraction: float = C.GIMME_BIOMASS_FRACTION,
              biomass_rxn: str = C.BIOMASS_RXN) -> GimmeResult:
    """Solve GIMME for one sample as a single weighted-parsimonious LP.

    Objective (minimise):  Σ_all (1 + c_i)·(fwd_i + rev_i)
    subject to steady state, the sample's exchange bounds, and biomass ≥
    fraction·biomass_max (RMF). The uniform weight 1 gives a unique parsimonious
    distribution (as in pFBA); the extra penalty c_i (>0 only for below-threshold
    reactions) pushes flux away from lowly-expressed reactions (GIMME). Reported
    ``inconsistency`` is the expression-penalty part Σ c_i·|v_i| at the optimum.
    Non-mutating (runs inside a ``with model`` context).
    """
    penalties, _ = _penalties(scores, quantile)

    with model:
        if exchange_bounds:
            _apply_exchange_bounds(model, exchange_bounds)

        # RMF: biomass must reach a fraction of its (constrained) optimum
        model.objective = biomass_rxn
        biomass_max = model.slim_optimize()
        if not np.isfinite(biomass_max) or biomass_max <= 0:
            return GimmeResult("infeasible_growth", float(biomass_max or 0.0),
                               np.nan, pd.Series(dtype=float))
        model.reactions.get_by_id(biomass_rxn).lower_bound = biomass_fraction * biomass_max

        # weighted-parsimonious objective over all reactions
        obj = model.problem.Objective(Zero, direction="min", name="gimme_pfba")
        model.objective = obj
        coeffs = {}
        for r in model.reactions:
            w = 1.0 + penalties.get(r.id, 0.0)
            coeffs[r.forward_variable] = w
            coeffs[r.reverse_variable] = w
        model.objective.set_linear_coefficients(coeffs)

        sol = model.optimize()
        if sol.status != "optimal":
            return GimmeResult(sol.status, np.nan, np.nan, pd.Series(dtype=float))

        fluxes = sol.fluxes
        inconsistency = float(sum(c * abs(fluxes.get(rid, 0.0))
                                  for rid, c in penalties.items()))
        growth = float(fluxes.get(biomass_rxn, np.nan))
        return GimmeResult("optimal", growth, inconsistency, fluxes, float(biomass_max))


def _apply_exchange_bounds(model, bounds: dict[str, tuple[float, float]]) -> None:
    ex_ids = {r.id for r in model.reactions}
    for rid, (lb, ub) in bounds.items():
        if rid in ex_ids:
            model.reactions.get_by_id(rid).bounds = (lb, ub)


def extract_panel(fluxes: pd.Series, model, subsystems=None, extra=None) -> pd.DataFrame:
    """Subset fluxes to the reporting panel, annotated with subsystem."""
    subsystems = set(subsystems or C.PANEL_SUBSYSTEMS)
    extra = set(extra or C.PANEL_EXTRA_RXNS)
    rows = []
    for rxn in model.reactions:
        sub = (rxn.subsystem or "").strip()
        if sub in subsystems or rxn.id in extra:
            rows.append({"reaction_id": rxn.id, "subsystem": sub or "OTHER",
                         "flux": float(fluxes.get(rxn.id, np.nan))})
    return pd.DataFrame(rows)
