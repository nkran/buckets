"""Assemble the long-format flux table + run manifest — the module deliverable."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C


def assemble_flux_table(records: list[dict]) -> pd.DataFrame:
    """``records`` = list of {project, clone, condition, timepoint, responder?,
    reaction_id, subsystem, flux}. Returns a tidy long DataFrame."""
    df = pd.DataFrame.from_records(records)
    key_order = [k for k in (C.SAMPLE_KEYS + ["responder"]) if k in df.columns]
    return df.sort_values(key_order + ["subsystem", "reaction_id"]).reset_index(drop=True)


def write_matrices(flux_cols: dict[str, pd.Series], expr_cols: dict[str, pd.Series],
                   out_dir: Path = C.FLUX_DIR) -> dict[str, Path]:
    """Write wide reaction × sample matrices for the Escher explorer.

    ``flux_cols`` / ``expr_cols`` map sample_id -> a per-reaction Series (full model
    fluxes and reaction expression scores respectively). Written as Parquet, with a
    CSV fallback if no parquet engine is available.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for name, cols in (("flux_matrix", flux_cols), ("rxn_expression", expr_cols)):
        if not cols:
            continue
        mat = pd.DataFrame(cols)
        mat.index.name = "reaction_id"
        try:
            p = out_dir / f"{name}.parquet"
            mat.to_parquet(p)
            paths[name] = p
        except Exception:
            p = out_dir / f"{name}.csv"
            mat.to_csv(p)
            paths[name] = p
    return paths


def write_outputs(flux_table: pd.DataFrame, manifest: dict,
                  out_dir: Path = C.FLUX_DIR) -> dict[str, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    csv_path = out_dir / "flux_table.csv"
    flux_table.to_csv(csv_path, index=False)
    paths["csv"] = csv_path
    try:
        pq_path = out_dir / "flux_table.parquet"
        flux_table.to_parquet(pq_path, index=False)
        paths["parquet"] = pq_path
    except Exception as exc:  # pyarrow missing / engine error — CSV still written
        manifest.setdefault("warnings", []).append(f"parquet export skipped: {exc}")
    man_path = out_dir / "manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2, default=_json_default))
    paths["manifest"] = man_path
    return paths


def build_manifest(*, model, solver: str, n_samples: int, coverage: dict,
                   infeasible: list[dict], params: dict) -> dict:
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "model": {
            "id": getattr(model, "id", None),
            "n_reactions": len(model.reactions),
            "n_metabolites": len(model.metabolites),
            "n_genes": len(model.genes),
            "biomass_reaction": C.BIOMASS_RXN,
        },
        "solver": solver,
        "n_samples": n_samples,
        "n_infeasible": len(infeasible),
        "infeasible_samples": infeasible,
        "mapping_coverage": coverage,
        "parameters": params,
    }


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, Path):
        return str(o)
    return str(o)
