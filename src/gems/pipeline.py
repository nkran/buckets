"""End-to-end orchestration: data → base model → q-rates → GIMME/FBA → export → viz.

Run:
    python -m gems.pipeline --synthetic            # full synthetic run
    python -m gems.pipeline --synthetic --limit 24 # quick subset
    python -m gems.pipeline --raw data/raw         # real data
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from . import config as C
from . import data_io, export, gene_mapping, gimme_fba, model_setup, qrate_calc, viz


# ---------------------------------------------------------------------------
# Per-sample work. Kept at module level so it can run in worker processes.
# ---------------------------------------------------------------------------
_W: dict = {}   # per-worker state: model + crosswalk, loaded once by _init_worker


def _init_worker(crosswalk: pd.DataFrame) -> None:
    _W["model"] = model_setup.get_base_model()
    _W["crosswalk"] = crosswalk


def _solve_sample(sid: str, expr: pd.DataFrame, rates: dict | None, g_meas: float,
                  qbounds: dict | None, model=None, crosswalk=None) -> dict:
    """Map expression, build exchange bounds, run GIMME for one sample.

    ``expr`` is the single-column counts table for this sample. Exactly one of
    ``rates`` (measured exchange rates, reconciled here) or ``qbounds`` (bounds from
    concentration-derived q-rates, or {} for expression-only) is used.
    """
    model = model if model is not None else _W["model"]
    crosswalk = crosswalk if crosswalk is not None else _W["crosswalk"]
    scores, rep = gene_mapping.map_expression(model, expr, crosswalk, sid)
    out = {"sid": sid, "scores": scores, "rep": rep, "recon": None}
    if rates is not None:
        bounds, recon = qrate_calc.reconcile_rates(model, rates, min_growth=g_meas)
        if bounds is None:
            out["status"] = "rates_not_reconcilable"
            return out
        out["recon"] = recon.assign(sample_id=sid)
    else:
        bounds = qbounds or {}
    res = gimme_fba.run_gimme(model, scores, exchange_bounds=bounds)
    out["status"] = res.status if not res.fluxes.empty else "no_fluxes"
    out["res"] = res
    if out["status"] == "optimal":
        out["panel"] = gimme_fba.extract_panel(res.fluxes, model)
    return out


def run(synthetic: bool = True, raw_dir=None, limit: int | None = None,
        rebuild_model: bool = False, make_figures: bool = True,
        contrasts=None, crosswalk_path=None, out_dir=None, synth_dir=None,
        jobs: int = 1) -> dict:
    """Run the whole pipeline. ``jobs`` > 1 solves samples in parallel worker
    processes (0 = all CPU cores); results are identical to a serial run."""
    # Synthetic runs never write into results/fluxes or results/figures (those hold the
    # real-data deliverable the explorer reads); they default to results/synthetic/.
    from pathlib import Path
    out_dir = Path(out_dir) if out_dir else (C.RESULTS / "synthetic" if synthetic else C.RESULTS)
    flux_dir, fig_dir = out_dir / "fluxes", out_dir / "figures"

    # 1. base CHO-S model
    model = model_setup.get_base_model(rebuild=rebuild_model)
    # optlang interface module name, e.g. "glpk_interface"
    solver = type(model.solver).__module__.rsplit(".", 1)[-1]

    # 2. data
    if synthetic:
        entrez = [g.id for g in model.genes]
        ds = data_io.make_synthetic_dataset(entrez, out_dir=synth_dir or C.SYNTH_DIR)
    else:
        ds = data_io.load_dataset(raw_dir or (C.DATA / "raw"))
    crosswalk = (ds.crosswalk if ds.crosswalk is not None
                 else gene_mapping.load_crosswalk(crosswalk_path or C.CROSSWALK_TSV))

    # 3. q-rates for all vessels (only if bioprocess measurements are present;
    #    transcriptomics-only datasets run GIMME with expression alone)
    have_bioprocess = all(c in ds.metadata.columns for c in C.BIOPROCESS_COLS)
    qrates = qrate_calc.compute_qrates(ds.metadata) if have_bioprocess else None
    # Alternative: measured rates supplied directly (mmol/gDW/h per exchange reaction)
    # in <raw_dir>/exchange_rates.csv, with optional <raw_dir>/measured_growth.csv.
    rates_table = growth_meas = None
    if not synthetic:
        rdir = Path(raw_dir or (C.DATA / "raw"))
        if (rdir / "exchange_rates.csv").exists():
            rates_table = pd.read_csv(rdir / "exchange_rates.csv", index_col=0)
            rates_table.index = rates_table.index.astype(str)
            if (rdir / "measured_growth.csv").exists():
                growth_meas = pd.read_csv(rdir / "measured_growth.csv", index_col=0).iloc[:, 0]
                growth_meas.index = growth_meas.index.astype(str)
            print(f"Using measured exchange rates for {rates_table.shape[1]} reactions")
    if not have_bioprocess and rates_table is None:
        print("No bioprocess columns found -> running expression-only (no exchange bounds)")
    recon_rows: list = []
    summary_rows: list[dict] = []

    # 4. per-sample GIMME + FBA
    sample_ids = ds.sample_ids()
    if limit:
        sample_ids = sample_ids[:limit]

    records: list[dict] = []
    infeasible: list[dict] = []
    cov_accum = {"entrez_coverage": [], "model_gene_coverage": [], "n_reactions_scored": []}
    flux_cols: dict[str, "pd.Series"] = {}   # sample_id -> full model fluxes (for explorer)
    expr_cols: dict[str, "pd.Series"] = {}   # sample_id -> reaction expression scores

    def task_args(sid):
        meta = ds.metadata.loc[sid]
        if rates_table is not None and sid in rates_table.index:
            g = float(growth_meas.get(sid, 0.0)) if growth_meas is not None else 0.0
            return (sid, ds.counts[[sid]], rates_table.loc[sid].to_dict(), g, None)
        qb = ({} if qrates is None else qrate_calc.qrate_bounds_for_sample(
            qrates, meta["project"], meta["clone"], meta["condition"], meta["timepoint"]))
        return (sid, ds.counts[[sid]], None, 0.0, qb)

    import os
    n_jobs = (os.cpu_count() or 1) if jobs == 0 else max(1, int(jobs))
    n_jobs = min(n_jobs, len(sample_ids)) or 1
    results: dict[str, dict] = {}
    if n_jobs == 1:
        for k, sid in enumerate(sample_ids, 1):
            results[sid] = _solve_sample(*task_args(sid), model=model, crosswalk=crosswalk)
            if k % 10 == 0 or k == len(sample_ids):
                print(f"[{k}/{len(sample_ids)}] {sid}: {results[sid]['status']}")
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        print(f"Solving {len(sample_ids)} samples on {n_jobs} worker processes")
        with ProcessPoolExecutor(max_workers=n_jobs, initializer=_init_worker,
                                 initargs=(crosswalk,)) as pool:
            futs = {pool.submit(_solve_sample, *task_args(sid)): sid for sid in sample_ids}
            for k, fut in enumerate(as_completed(futs), 1):
                r = fut.result()
                results[r["sid"]] = r
                if k % 10 == 0 or k == len(sample_ids):
                    print(f"[{k}/{len(sample_ids)}] done")

    # collect in sample order so output is identical to a serial run
    for sid in sample_ids:
        r = results[sid]
        meta = ds.metadata.loc[sid]
        rep = r["rep"]
        cov_accum["entrez_coverage"].append(rep.entrez_coverage)
        cov_accum["model_gene_coverage"].append(rep.model_gene_coverage)
        cov_accum["n_reactions_scored"].append(rep.n_reactions_scored)
        if r["status"] != "optimal":
            infeasible.append({"sample_id": sid, "status": r["status"]})
            print(f"{sid}: {r['status']}")
            continue
        res = r["res"]
        if r["recon"] is not None:
            recon_rows.append(r["recon"])
        summary_rows.append({
            "sample_id": sid, "project": meta["project"], "clone": meta["clone"],
            "condition": meta["condition"], "timepoint": meta["timepoint"],
            "growth_capacity": res.growth_max, "growth_used": res.growth,
            "measured_growth": (float(growth_meas.get(sid, np.nan))
                                if growth_meas is not None else np.nan),
            "expression_inconsistency": res.inconsistency,
            "rate_adjustment_total": (float(r["recon"]["deviation"].abs().sum())
                                      if r["recon"] is not None else np.nan),
        })
        # full flux + expression vectors for the interactive explorer
        flux_cols[sid] = res.fluxes
        expr_cols[sid] = pd.Series(r["scores"], dtype=float)
        for _, prow in r["panel"].iterrows():
            records.append({
                "project": meta["project"], "clone": meta["clone"],
                "condition": meta["condition"], "timepoint": meta["timepoint"],
                "responder": meta.get("responder", np.nan),
                "reaction_id": prow["reaction_id"], "subsystem": prow["subsystem"],
                "flux": prow["flux"],
            })

    # 5. export
    flux_table = export.assemble_flux_table(records)
    coverage = {kk: (round(float(np.mean(vv)), 4) if vv else None)
                for kk, vv in cov_accum.items()}
    manifest = export.build_manifest(
        model=model, solver=solver, n_samples=len(sample_ids), coverage=coverage,
        infeasible=infeasible,
        params={"gimme_quantile": C.GIMME_QUANTILE,
                "gimme_biomass_fraction": C.GIMME_BIOMASS_FRACTION,
                "qrate_tol": qrate_calc.QRATE_TOL,
                "exchange_constraints_applied": bool(have_bioprocess or rates_table is not None),
                "exchange_constraint_source": ("measured rates (reconciled)" if rates_table is not None
                                               else "q-rates from concentrations" if have_bioprocess
                                               else "none"),
                "panel_subsystems": C.PANEL_SUBSYSTEMS},
    )
    paths = export.write_outputs(flux_table, manifest, out_dir=flux_dir)
    mat_paths = export.write_matrices(flux_cols, expr_cols, out_dir=flux_dir)
    pd.DataFrame(summary_rows).to_csv(flux_dir / "sample_summary.csv", index=False)
    if recon_rows:
        pd.concat(recon_rows)[["sample_id", "reaction_id", "measured", "used", "deviation"]] \
            .to_csv(flux_dir / "exchange_reconciliation.csv", index=False)
    print(f"\nWrote flux table ({flux_table.shape[0]} rows) -> {paths.get('csv')}")
    if mat_paths:
        print(f"Wrote explorer matrices -> {', '.join(p.name for p in mat_paths.values())}")
    print(f"Coverage: {coverage} | infeasible: {len(infeasible)}/{len(sample_ids)}")

    # 6. figures
    figs = []
    if make_figures and not flux_table.empty:
        figs = viz.make_all_figures(flux_table, out_dir=fig_dir, contrasts=contrasts)
        print(f"Wrote {len(figs)} figures -> {fig_dir}")

    return {"flux_table": flux_table, "manifest": manifest, "paths": paths, "figures": figs}


def main():
    ap = argparse.ArgumentParser(description="CHO iFB mechanistic flux-generation pipeline")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--synthetic", action="store_true", help="use generated placeholder data")
    g.add_argument("--raw", metavar="DIR", help="directory with counts.csv + metadata.csv")
    ap.add_argument("--limit", type=int, default=None, help="cap number of samples")
    ap.add_argument("--rebuild-model", action="store_true", help="rebuild cached base model")
    ap.add_argument("--crosswalk", metavar="TSV", default=None,
                    help="gene crosswalk tsv (ensembl_gene_id, entrez_id); "
                         "default data/maps/ensembl_picrh_to_entrez.tsv")
    ap.add_argument("--out", metavar="DIR", default=None,
                    help="output root (default: results/ for real data, results/synthetic/ for --synthetic)")
    ap.add_argument("--jobs", type=int, default=0,
                    help="parallel worker processes (default 0 = all CPU cores; 1 = serial)")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args()

    synthetic = args.synthetic or not args.raw
    run(synthetic=synthetic, raw_dir=args.raw, limit=args.limit,
        rebuild_model=args.rebuild_model, make_figures=not args.no_figures,
        crosswalk_path=args.crosswalk, out_dir=args.out, jobs=args.jobs)


if __name__ == "__main__":
    main()
