"""Prepare the ZeLa fed-batch example (the dataset published with iCHO3K).

Study: wild-type CHO vs "ZeLa" (a CHO line engineered to secrete no lactate), grown in
fed-batch bioreactors and sampled on days 2, 4, 6, 8, 12 and 14. Source files live in
https://github.com/LewisLabUCSD/iCHO3K (Data/ZeLa Data, Data/Uptake_Secretion_Rates).

Unlike GSE140669 this dataset has BOTH layers the model can use:
  * transcriptomics  - TPM per gene (symbols) for 43 samples
  * exchange rates   - measured uptake/secretion per bioreactor and interval, already in
                       mmol/gDW/h and already named by iCHO3K exchange reaction id

Writes into data/zela/:
  counts.csv          genes x samples (TPM)
  metadata.csv        sample -> project, clone (bioreactor), condition (WT/ZeLa), timepoint (day)
  exchange_rates.csv  samples x exchange reactions (mmol/gDW/h; negative = uptake)
  measured_growth.csv sample -> measured specific growth rate (1/h), kept for comparison only
  crosswalk.tsv       gene symbol -> iCHO3K Entrez id

Run:   python scripts/prepare_zela.py
Then:  python -m gems.pipeline --raw data/zela --crosswalk data/zela/crosswalk.tsv --out results/zela
"""

from __future__ import annotations

import base64
import json
import urllib.request

import pandas as pd

from gems import config as C
from gems import gene_ids as gid

OUT = C.DATA / "zela"
REPO_API = "https://api.github.com/repos/LewisLabUCSD/iCHO3K"
FILES = {
    "Data/ZeLa Data/zela_transcriptomics_data_tpm.tsv": "tpm.tsv",
    "Data/Uptake_Secretion_Rates/WT_pkl_mat_excel/uptake_secretion_rates_wt.csv": "rates_wt.csv",
    "Data/Uptake_Secretion_Rates/ZeLa_pkl_mat_excel/uptake_secretion_rates_zela.csv": "rates_zela.csv",
}

# Transcribed from "Data/ZeLa Data/Bio141-148_Sample_info.numbers" in the iCHO3K repo:
# bioreactor -> (strain, first sample number, days sampled). Sample numbers run
# consecutively (S1..S43) in this order.
BIOREACTORS = [
    ("U1", "WT", [2, 4, 6, 8, 12, 14]),
    ("U2", "WT", [2, 4, 6, 8, 12, 14]),
    ("U3", "WT", [2, 4, 6, 8, 14]),
    ("U4", "ZeLa", [4, 6, 8, 14]),
    ("U5", "ZeLa", [2, 4, 8, 12, 14]),
    ("U6", "ZeLa", [2, 4, 6, 8, 12, 14]),
    ("U7", "ZeLa", [2, 4, 6, 8, 14]),
    ("U8", "ZeLa", [2, 4, 6, 8, 12, 14]),
]
# A sample taken on day D is matched to the rate interval that ENDS on day D.
INTERVAL_END_DAY = {"P0 to P2": 2, "P2 to P4": 4, "P4 to P6": 6, "P6 to P8": 8,
                    "P8 to P12": 12, "P12 to P14": 14}


def download() -> None:
    """Fetch the source files through the git-blob API (raw.githubusercontent rate-limits)."""
    OUT.mkdir(parents=True, exist_ok=True)
    if all((OUT / name).exists() for name in FILES.values()):
        return
    tree = json.load(urllib.request.urlopen(f"{REPO_API}/git/trees/main?recursive=1"))
    for item in tree["tree"]:
        name = FILES.get(item["path"])
        if name and not (OUT / name).exists():
            blob = json.load(urllib.request.urlopen(f"{REPO_API}/git/blobs/{item['sha']}"))
            (OUT / name).write_bytes(base64.b64decode(blob["content"]))
            print(f"downloaded {name}")


def sample_table() -> pd.DataFrame:
    rows, n = [], 0
    for reactor, strain, days in BIOREACTORS:
        for day in days:
            n += 1
            rows.append({"sample_id": f"S{n}", "project": "ZeLa", "clone": reactor,
                         "condition": strain, "timepoint": day})
    return pd.DataFrame(rows).set_index("sample_id")


def main() -> None:
    download()
    meta = sample_table()

    # 1. expression: TPM, genes x samples, columns in S1..S43 order
    tpm = pd.read_csv(OUT / "tpm.tsv", sep="\t", index_col=0)
    tpm.index.name = "gene"
    tpm = tpm[list(meta.index)]
    tpm.to_csv(OUT / "counts.csv")

    # 2. measured exchange rates: one row per sample (bioreactor + interval ending that day)
    rates = pd.concat([pd.read_csv(OUT / "rates_wt.csv"), pd.read_csv(OUT / "rates_zela.csv")])
    rates["day"] = rates["Interval"].map(INTERVAL_END_DAY)
    rates = rates.set_index(["Batch ID", "day"])
    ex_cols = [c for c in rates.columns if c.startswith("EX_")]
    key = list(zip(meta["clone"], meta["timepoint"]))
    per_sample = rates.loc[key]
    ex = per_sample[ex_cols].copy()
    ex.index = meta.index
    ex.to_csv(OUT / "exchange_rates.csv")
    growth = per_sample["exp_growth_rate"].copy()
    growth.index = meta.index
    growth.rename("measured_growth_per_h").to_csv(OUT / "measured_growth.csv")

    meta.to_csv(OUT / "metadata.csv")

    # 3. gene crosswalk (symbols -> iCHO3K Entrez) from the master gene-id bridge
    gid.make_best_crosswalk(tpm.index, OUT / "crosswalk.tsv")

    print(f"samples: {len(meta)}  ({(meta.condition == 'WT').sum()} WT, "
          f"{(meta.condition == 'ZeLa').sum()} ZeLa)")
    print(f"genes: {tpm.shape[0]}   exchange reactions with measured rates: {len(ex_cols)}")
    print(f"wrote counts.csv, metadata.csv, exchange_rates.csv, measured_growth.csv, "
          f"crosswalk.tsv -> {OUT}")


if __name__ == "__main__":
    main()
