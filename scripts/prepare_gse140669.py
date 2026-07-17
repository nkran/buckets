"""Prepare a real public CHO RNA-seq dataset (GEO GSE140669) for the gems pipeline.

Study: "Defining lncRNAs Associated with CHO Cell Growth and IgG Productivity by
RNA-Seq" (GSE140669). 4 CHO clones (3936, 3068, 4384, 3080) × 3 timepoints
(ST = seed/standard, D4 = day 4, D7 = day 7 of culture) × 3 technical replicates.

This script:
  1. collapses the technical triplicates per (clone, timepoint) by summing counts
     (as the study design note specifies: triplicates -> one high-depth reading),
  2. CPM-normalizes each collapsed library (stand-in for DESeq2 normalization so
     sequencing-depth differences don't bias GIMME),
  3. writes data/raw/counts.csv + data/raw/metadata.csv in the pipeline schema,
  4. builds the gene crosswalk (dataset symbol/LOC -> iCHO3K Entrez) from the
     master gene-id bridge, written to data/maps/gse140669_symbol_to_entrez.tsv.

Run:  python scripts/prepare_gse140669.py
Then: python -m gems.pipeline --raw data/raw
"""

from __future__ import annotations

import pandas as pd

from gems import config as C

RAW_COUNTS = C.DATA / "raw" / "GSE140669_counts.csv"
MODEL_MAP = C.DATA / "maps" / "icho3k_gene_id_mapping.csv"

TIMEPOINT_DAY = {"ST": 0, "D4": 4, "D7": 7}   # culture day


def main() -> None:
    counts = pd.read_csv(RAW_COUNTS, index_col=0)
    counts.index.name = "gene"

    # --- 1. collapse technical replicates: sum over (clone, timepoint) ---
    triples = {}          # (clone, tp) -> [sample cols]
    for col in counts.columns:
        clone, tp, _rep = col.split("_")
        triples.setdefault((clone, tp), []).append(col)

    collapsed = {}
    for (clone, tp), cols in triples.items():
        collapsed[f"{clone}_{tp}"] = counts[cols].sum(axis=1)
    mat = pd.DataFrame(collapsed)

    # --- 2. CPM normalization (per collapsed library) ---
    cpm = mat.div(mat.sum(axis=0), axis=1) * 1e6

    # --- 3. metadata ---
    rows = []
    for sid in cpm.columns:
        clone, tp = sid.split("_")
        rows.append({"sample_id": sid, "project": "GSE140669", "clone": clone,
                     "condition": "FB", "timepoint": TIMEPOINT_DAY[tp]})
    metadata = pd.DataFrame(rows).set_index("sample_id")

    (C.DATA / "raw").mkdir(parents=True, exist_ok=True)
    cpm.to_csv(C.DATA / "raw" / "counts.csv")
    metadata.to_csv(C.DATA / "raw" / "metadata.csv")

    # --- 4. gene crosswalk: dataset gene -> iCHO3K Entrez ---
    # GSE140669 is keyed by NCBI symbol/LOC (not Ensembl PICRH). Build the best crosswalk
    # from the master gene-id bridge (symbol -> Entrez), unioned with iCHO3K's own symbol
    # map. This reaches ~89% model-gene coverage vs ~66% from the iCHO3K map alone.
    from gems import gene_ids as gid
    xwalk_path = C.DATA / "maps" / "gse140669_symbol_to_entrez.tsv"
    if (C.DATA / "maps" / "gene_id_master.tsv").exists():
        gid.make_best_crosswalk(cpm.index, xwalk_path)
    else:   # fallback if the master hasn't been built yet
        m = pd.read_csv(MODEL_MAP, dtype=str)
        m["entrez_gene_id"] = m["entrez_gene_id"].str.split(".").str[0]
        (m.rename(columns={"gene": "ensembl_gene_id", "entrez_gene_id": "entrez_id"})
          [["ensembl_gene_id", "entrez_id"]].dropna().drop_duplicates()
          .to_csv(xwalk_path, sep="\t", index=False))
    xwalk = pd.read_csv(xwalk_path, sep="\t", dtype=str)

    matched = counts.index.intersection(xwalk["ensembl_gene_id"])
    print(f"collapsed samples : {list(cpm.columns)}")
    print(f"genes             : {cpm.shape[0]}  (matched to crosswalk: {len(matched)})")
    print(f"wrote counts.csv, metadata.csv -> {C.DATA/'raw'}")
    print(f"wrote crosswalk ({len(xwalk)} rows) -> {xwalk_path}")
    print(f"\nRun with: python -m gems.pipeline --raw data/raw --crosswalk {xwalk_path}")


if __name__ == "__main__":
    main()
