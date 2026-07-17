"""Load / validate transcriptomic + bioprocess inputs, and generate a schema-correct
synthetic dataset so the pipeline runs end-to-end before real data is available.

Real inputs expected (drop into ``data/raw/``):
  * ``counts.csv``    — DESeq2 *normalized* counts, genes (Ensembl PICRH ``ENSCGRG…``)
                        as rows, sample ids as columns.
  * ``metadata.csv``  — one row per sample id, columns:
                        project, clone, condition, timepoint + bioprocess measurements
                        (see ``config.BIOPROCESS_COLS``). Optional ``responder`` label.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C


@dataclass
class Dataset:
    counts: pd.DataFrame       # genes (Ensembl PICRH) x samples, normalized
    metadata: pd.DataFrame     # index = sample_id, columns incl. SAMPLE_KEYS + BIOPROCESS_COLS
    crosswalk: pd.DataFrame | None = None  # ensembl_gene_id, entrez_id (synthetic runs)

    def sample_ids(self) -> list[str]:
        return list(self.metadata.index)


def sample_id(project, clone, condition, timepoint) -> str:
    return f"{project}_C{clone:02d}_{condition}_D{timepoint:02d}"


# ---------------------------------------------------------------------------
# Loading + validation
# ---------------------------------------------------------------------------
def load_dataset(raw_dir: Path = C.DATA / "raw") -> Dataset:
    raw_dir = Path(raw_dir)
    counts = pd.read_csv(raw_dir / "counts.csv", index_col=0)
    metadata = pd.read_csv(raw_dir / "metadata.csv", index_col=0)
    validate(counts, metadata)
    return Dataset(counts=counts, metadata=metadata)


def validate(counts: pd.DataFrame, metadata: pd.DataFrame) -> None:
    missing_keys = [k for k in C.SAMPLE_KEYS if k not in metadata.columns]
    if missing_keys:
        raise ValueError(f"metadata missing required columns: {missing_keys}")
    if not set(counts.columns) >= set(metadata.index):
        orphan = set(metadata.index) - set(counts.columns)
        raise ValueError(f"{len(orphan)} samples in metadata have no counts column: "
                         f"{sorted(orphan)[:5]}…")
    if (counts.values < 0).any():
        raise ValueError("counts contain negative values (expected DESeq2 normalized counts)")
    if counts.index.duplicated().any():
        raise ValueError("duplicate gene ids in counts index")


# ---------------------------------------------------------------------------
# Synthetic dataset
# ---------------------------------------------------------------------------
def make_synthetic_dataset(
    entrez_genes: list[str],
    out_dir: Path = C.SYNTH_DIR,
    n_extra_unmapped: int = 300,
    seed: int = 0,
) -> Dataset:
    """Build a placeholder dataset keyed by fake Ensembl PICRH gene ids.

    ``entrez_genes`` are the model's Entrez ids; each is given a synthetic
    ``ENSCGRG…`` accession so gene mapping can be exercised. A fraction of model
    genes is deliberately left unmapped, and extra unmapped Ensembl genes are
    added, to mimic realistic crosswalk coverage loss.
    """
    rng = np.random.default_rng(seed)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- crosswalk: map ~85% of model genes to a synthetic Ensembl id ---
    entrez_genes = list(entrez_genes)
    n_mapped = int(0.85 * len(entrez_genes))
    mapped_entrez = list(rng.permutation(entrez_genes)[:n_mapped])
    ens_ids = [f"ENSCGRG00000{ i:06d}" for i in range(1, n_mapped + 1)]
    crosswalk = pd.DataFrame({"ensembl_gene_id": ens_ids, "entrez_id": mapped_entrez})

    # extra expressed genes with no model counterpart (coverage-loss realism)
    extra_ens = [f"ENSCGRG00000{ i:06d}" for i in range(n_mapped + 1, n_mapped + 1 + n_extra_unmapped)]
    all_ens = ens_ids + extra_ens

    # --- design / metadata ---
    rows = []
    # assign each clone a latent "responder" status (per project, per clone)
    responder = {}
    for p in range(1, C.N_PROJECTS + 1):
        for c in range(1, C.N_CLONES + 1):
            responder[(p, c)] = bool(rng.random() < 0.5)

    for p in range(1, C.N_PROJECTS + 1):
        proj = f"P{p}"
        for c in range(1, C.N_CLONES + 1):
            is_resp = responder[(p, c)]
            for cond in C.CONDITIONS:
                traj = _bioprocess_trajectory(cond, is_resp, rng)
                for i, day in enumerate(C.TIMEPOINTS):
                    sid = sample_id(proj, c, cond, day)
                    rows.append({
                        "sample_id": sid, "project": proj, "clone": c,
                        "condition": cond, "timepoint": day,
                        "responder": "Responder" if is_resp else "NonResponder",
                        **{col: traj[col][i] for col in C.BIOPROCESS_COLS},
                    })
    metadata = pd.DataFrame(rows).set_index("sample_id")

    # --- expression: structured lognormal with condition/phenotype signal ---
    counts = _synthetic_expression(all_ens, metadata, rng)

    # persist
    counts.to_csv(out_dir / "counts.csv")
    metadata.to_csv(out_dir / "metadata.csv")
    # kept beside the synthetic data: must never overwrite the real BioMart crosswalk
    crosswalk.to_csv(out_dir / "crosswalk.tsv", sep="\t", index=False)

    validate(counts, metadata)
    return Dataset(counts=counts, metadata=metadata, crosswalk=crosswalk)


def _bioprocess_trajectory(condition: str, responder: bool, rng) -> dict:
    """Return per-timepoint arrays for each bioprocess column (length = TIMEPOINTS)."""
    n = len(C.TIMEPOINTS)
    t = np.array(C.TIMEPOINTS, dtype=float)
    intensified = condition == "iFB"

    # VCD: intensified starts higher (high inoculation), grows then declines.
    vcd0 = 8.0 if intensified else 3.0
    peak_day = 9.0 if responder else 7.0
    vcd = vcd0 + (14.0 if responder else 9.0) * np.exp(-((t - peak_day) ** 2) / (2 * 3.0 ** 2))
    # Non-responders under iFB crash late
    if intensified and not responder:
        vcd = vcd * np.clip(1.0 - 0.10 * (t - peak_day), 0.2, 1.0)
    vcd *= 1 + 0.03 * rng.standard_normal(n)

    # Viability declines; non-responder iFB crashes.
    viab = np.clip(99 - 0.6 * (t - 3), 60, 100)
    if intensified and not responder:
        viab = np.clip(99 - 2.8 * (t - 3), 30, 100)
    viab *= 1 + 0.01 * rng.standard_normal(n)

    # Glucose consumed over time (intensified consumes faster).
    g0 = 45.0
    gluc = np.clip(g0 - (3.2 if intensified else 2.2) * (t - 3), 2, g0)
    gluc *= 1 + 0.02 * rng.standard_normal(n)

    # Glutamine consumed.
    gln = np.clip(6.0 - 0.45 * (t - 3), 0.2, 6.0) * (1 + 0.02 * rng.standard_normal(n))
    # Glutamate mild accumulation.
    glu = (1.0 + 0.15 * (t - 3)) * (1 + 0.03 * rng.standard_normal(n))

    # Lactate: rises early; responders shift to CONSUMPTION (metabolic switch).
    lac = 2.0 + 1.6 * (t - 3)
    if responder:
        # peak then decline (lactate consumption)
        lac = 2.0 + 6.0 * np.exp(-((t - 7) ** 2) / (2 * 3.0 ** 2)) + 0.4 * (t - 3)
    if intensified and not responder:
        lac = 2.0 + 2.4 * (t - 3)  # keeps climbing
    lac = np.clip(lac, 0.5, None) * (1 + 0.03 * rng.standard_normal(n))

    # Ammonia accumulates.
    nh3 = (0.5 + 0.35 * (t - 3)) * (1 + 0.03 * rng.standard_normal(n))

    # Titer accumulates; responders titre higher, esp. under iFB.
    rate = 90 if responder else 55
    if intensified:
        rate *= 1.25 if responder else 0.9
    titer = np.clip(rate * (t - 3), 0, None) * (1 + 0.03 * rng.standard_normal(n))

    return {"VCD": vcd, "Viability": viab, "Titer": titer, "Gluc": gluc,
            "Gln": gln, "Glu": glu, "Lac": lac, "NH3": nh3}


def _synthetic_expression(gene_ids: list[str], metadata: pd.DataFrame, rng) -> pd.DataFrame:
    """Lognormal baseline with a UPR/stress module up in iFB non-responders and a
    metabolic-robustness module up in responders — gives the flux layer real signal."""
    n_genes = len(gene_ids)
    samples = list(metadata.index)
    # baseline per-gene mean expression (normalized counts scale)
    base_mu = rng.lognormal(mean=4.0, sigma=1.2, size=n_genes)
    mat = np.outer(base_mu, np.ones(len(samples)))
    mat *= rng.lognormal(mean=0.0, sigma=0.25, size=mat.shape)  # sample noise

    # designate gene modules
    idx = rng.permutation(n_genes)
    upr = idx[: n_genes // 20]                 # ~5% stress/UPR
    robust = idx[n_genes // 20: n_genes // 10]  # next ~5% robustness/oxphos

    for j, sid in enumerate(samples):
        row = metadata.loc[sid]
        intensified = row["condition"] == "iFB"
        is_resp = row["responder"] == "Responder"
        if intensified and not is_resp:
            mat[upr, j] *= 2.5
        if is_resp:
            mat[robust, j] *= 1.8
        # late-timepoint amplification of the divergence
        if row["timepoint"] >= 9:
            if intensified and not is_resp:
                mat[upr, j] *= 1.4
            if is_resp:
                mat[robust, j] *= 1.2

    return pd.DataFrame(mat, index=gene_ids, columns=samples)
