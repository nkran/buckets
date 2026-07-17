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


def run(synthetic: bool = True, raw_dir=None, limit: int | None = None,
        rebuild_model: bool = False, make_figures: bool = True,
        contrasts=None, crosswalk_path=None, out_dir=None, synth_dir=None) -> dict:
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
    if not have_bioprocess:
        print("No bioprocess columns found -> running expression-only (no exchange bounds)")

    # 4. per-sample GIMME + FBA
    sample_ids = ds.sample_ids()
    if limit:
        sample_ids = sample_ids[:limit]

    records: list[dict] = []
    infeasible: list[dict] = []
    cov_accum = {"entrez_coverage": [], "model_gene_coverage": [], "n_reactions_scored": []}
    flux_cols: dict[str, "pd.Series"] = {}   # sample_id -> full model fluxes (for explorer)
    expr_cols: dict[str, "pd.Series"] = {}   # sample_id -> reaction expression scores

    for k, sid in enumerate(sample_ids, 1):
        meta = ds.metadata.loc[sid]
        scores, rep = gene_mapping.map_expression(model, ds.counts, crosswalk, sid)
        cov_accum["entrez_coverage"].append(rep.entrez_coverage)
        cov_accum["model_gene_coverage"].append(rep.model_gene_coverage)
        cov_accum["n_reactions_scored"].append(rep.n_reactions_scored)

        bounds = ({} if qrates is None else qrate_calc.qrate_bounds_for_sample(
            qrates, meta["project"], meta["clone"], meta["condition"], meta["timepoint"]))
        res = gimme_fba.run_gimme(model, scores, exchange_bounds=bounds)

        if res.status != "optimal" or res.fluxes.empty:
            infeasible.append({"sample_id": sid, "status": res.status})
            print(f"[{k}/{len(sample_ids)}] {sid}: INFEASIBLE ({res.status})")
            continue

        # full flux + expression vectors for the interactive explorer
        flux_cols[sid] = res.fluxes
        expr_cols[sid] = pd.Series(scores, dtype=float)

        panel = gimme_fba.extract_panel(res.fluxes, model)
        for _, prow in panel.iterrows():
            records.append({
                "project": meta["project"], "clone": meta["clone"],
                "condition": meta["condition"], "timepoint": meta["timepoint"],
                "responder": meta.get("responder", np.nan),
                "reaction_id": prow["reaction_id"], "subsystem": prow["subsystem"],
                "flux": prow["flux"],
            })
        if k % 10 == 0 or k == len(sample_ids):
            print(f"[{k}/{len(sample_ids)}] {sid}: growth={res.growth:.3f} "
                  f"incon={res.inconsistency:.2f}")

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
                "exchange_constraints_applied": bool(have_bioprocess),
                "panel_subsystems": C.PANEL_SUBSYSTEMS},
    )
    paths = export.write_outputs(flux_table, manifest, out_dir=flux_dir)
    mat_paths = export.write_matrices(flux_cols, expr_cols, out_dir=flux_dir)
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
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args()

    synthetic = args.synthetic or not args.raw
    run(synthetic=synthetic, raw_dir=args.raw, limit=args.limit,
        rebuild_model=args.rebuild_model, make_figures=not args.no_figures,
        crosswalk_path=args.crosswalk, out_dir=args.out)


if __name__ == "__main__":
    main()
