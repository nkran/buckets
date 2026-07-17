"""Load iCHO3K and prepare the base CHO-S cell-line model.

The base model = iCHO3K (unblocked) with the CHO-S biomass objective and a defined
culture medium opened for uptake, verified to grow. It is the reusable description
of our sampled cell line; per-sample context models are derived from it.
"""

from __future__ import annotations

from pathlib import Path

import cobra

from . import config as C


def load_raw_model(path: Path = C.MODEL_JSON) -> cobra.Model:
    model = cobra.io.load_json_model(str(path))
    # open-source LP solver; GIMME/FBA are LP so this is sufficient.
    for solver in ("glpk", "hybrid"):
        try:
            model.solver = solver
            break
        except Exception:
            continue
    return model


def apply_defined_medium(model: cobra.Model) -> cobra.Model:
    """Close all exchanges, then open the defined base medium for uptake.

    Modifies ``model`` in place and returns it.
    """
    ex_ids = {r.id for r in model.exchanges}
    for ex in model.exchanges:
        # keep secretion open (ub>=0), close uptake by default
        ex.lower_bound = 0.0
    for rid in C.BASE_MEDIUM_UPTAKES:
        if rid in ex_ids:
            model.reactions.get_by_id(rid).lower_bound = -C.MEDIUM_UPTAKE
    for rid in C.MEDIUM_FREE:
        if rid in ex_ids:
            model.reactions.get_by_id(rid).lower_bound = -1000.0
    return model


def build_base_model(
    raw_path: Path = C.MODEL_JSON,
    out_path: Path = C.BASE_MODEL_JSON,
    cache: bool = True,
) -> cobra.Model:
    """Return the curated CHO-S base model, building + caching it if needed."""
    model = load_raw_model(raw_path)
    if C.BIOMASS_RXN not in {r.id for r in model.reactions}:
        raise ValueError(f"biomass reaction {C.BIOMASS_RXN!r} not found in model")
    model.objective = C.BIOMASS_RXN
    apply_defined_medium(model)

    sol = model.optimize()
    if sol.status != "optimal" or (sol.objective_value or 0) <= C.CHO_S_GROWTH_TOL:
        raise RuntimeError(
            f"Base CHO-S model does not grow on the defined medium "
            f"(status={sol.status}, growth={sol.objective_value}). "
            f"Revisit BASE_MEDIUM_UPTAKES in config.")
    model.id = "iCHO3K_base_cho_s"
    if cache:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        cobra.io.save_json_model(model, str(out_path))
    return model


def get_base_model(rebuild: bool = False) -> cobra.Model:
    """Load the cached base model, building it on first use."""
    if C.BASE_MODEL_JSON.exists() and not rebuild:
        model = cobra.io.load_json_model(str(C.BASE_MODEL_JSON))
        try:
            model.solver = "glpk"
        except Exception:
            pass
        return model
    return build_base_model()


def baseline_growth(model: cobra.Model) -> float:
    return float(model.slim_optimize())
