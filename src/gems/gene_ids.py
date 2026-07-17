"""Accessor for the master gene-identifier bridge (`data/maps/gene_id_master.tsv`).

The master table links, per Chinese-hamster (CriGri-PICRH-1.0) gene, the namespaces
that appear across CHO transcriptomics + metabolic modelling:

    ensembl_gene_id · entrez_id · symbol · refseq_mrna · refseq_protein ·
    ensembl_protein · uniprot · gene_biotype · chromosome · in_icho3k · icho3k_n_reactions

Multi-valued columns (`refseq_*`, `ensembl_protein`, `uniprot`) are ``;``-joined lists.

Typical uses::

    from gems import gene_ids as gid
    m = gid.load_master()                      # the whole table
    gid.translate(['Has2'], 'symbol', 'entrez_id')          # {'Has2': ['100751055']}
    gid.make_pipeline_crosswalk('refseq_mrna', 'my_xwalk.tsv')  # for a RefSeq-keyed dataset
    gid.detect_namespace(counts.index)          # guess how a matrix is keyed
"""

from __future__ import annotations

import functools
from pathlib import Path

import pandas as pd

from . import config as C

MASTER_TSV = C.DATA / "maps" / "gene_id_master.tsv"

# columns that hold ``;``-separated lists of ids
MULTI_COLS = {"refseq_mrna", "refseq_protein", "ensembl_protein", "uniprot"}
NAMESPACES = ["ensembl_gene_id", "entrez_id", "symbol", "refseq_mrna",
              "refseq_protein", "ensembl_protein", "uniprot"]


@functools.lru_cache(maxsize=1)
def load_master(path: str | Path = MASTER_TSV) -> pd.DataFrame:
    """Load the master bridge table (cached)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    if "in_icho3k" in df.columns:
        df["in_icho3k"] = df["in_icho3k"].map({"True": True, "False": False}).fillna(False)
    return df


def _exploded(source: str, target: str) -> pd.DataFrame:
    """Long (source_id, target_id) pairs, expanding any ``;``-joined list columns."""
    m = load_master()
    if source == target:                       # identity mapping
        vals = m[source].dropna()
        if source in MULTI_COLS:
            vals = vals.str.split(";").explode()
        vals = vals[vals.astype(str).str.len() > 0].drop_duplicates()
        return pd.DataFrame({source: vals.values, target: vals.values})
    cols = m[[source, target]].dropna()
    for col in (source, target):
        if col in MULTI_COLS:
            cols = cols.assign(**{col: cols[col].str.split(";")}).explode(col)
    cols = cols[(cols[source].astype(str).str.len() > 0)
                & (cols[target].astype(str).str.len() > 0)]
    return cols.drop_duplicates()


def crosswalk(source: str, target: str = "entrez_id") -> pd.DataFrame:
    """Two-column DataFrame mapping `source` namespace -> `target` namespace."""
    if source not in NAMESPACES or target not in NAMESPACES:
        raise ValueError(f"namespace must be one of {NAMESPACES}")
    return _exploded(source, target).reset_index(drop=True)


def translate(ids, source: str, target: str = "entrez_id") -> dict[str, list[str]]:
    """Map an iterable of ids from `source` to a list of `target` ids each."""
    x = crosswalk(source, target)
    lut: dict[str, list[str]] = {}
    for s, t in zip(x[source], x[target]):
        lut.setdefault(s, []).append(t)
    return {i: lut.get(str(i), []) for i in ids}


def make_best_crosswalk(ids, out_path: str | Path, *, include_icho3k_map: bool = True) -> Path:
    """Auto-detect a dataset's gene namespace and write the best entrez crosswalk.

    Uses the detected namespace's master mapping; for symbol-keyed data it also unions
    iCHO3K's own `gene_id_mapping.csv` (another symbol source) for maximum recall.
    Output columns are `ensembl_gene_id, entrez_id` (pipeline-compatible; the first
    column holds the source ids whatever their namespace).
    """
    ns = detect_namespace(ids)
    if ns is None:
        raise ValueError("could not detect gene id namespace from the given ids")
    x = crosswalk(ns, "entrez_id").rename(columns={ns: "ensembl_gene_id"})
    if ns == "symbol" and include_icho3k_map:
        icho = C.DATA / "maps" / "icho3k_gene_id_mapping.csv"
        if icho.exists():
            extra = pd.read_csv(icho, dtype=str)
            extra["entrez_gene_id"] = extra["entrez_gene_id"].str.split(".").str[0]
            extra = extra.rename(columns={"gene": "ensembl_gene_id",
                                          "entrez_gene_id": "entrez_id"})[["ensembl_gene_id", "entrez_id"]]
            x = pd.concat([x, extra]).dropna().drop_duplicates()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    x.to_csv(out_path, sep="\t", index=False)
    return out_path


def make_pipeline_crosswalk(source: str, out_path: str | Path,
                            target: str = "entrez_id") -> Path:
    """Write a `gems.pipeline`-compatible crosswalk (cols: ensembl_gene_id, entrez_id).

    The first column is named ``ensembl_gene_id`` for pipeline compatibility but holds
    ids of whatever `source` namespace you choose — so a dataset keyed by RefSeq or gene
    symbol can be mapped through the same code path. Pass the file via
    ``python -m gems.pipeline --raw … --crosswalk <out_path>``.
    """
    x = crosswalk(source, target).rename(
        columns={source: "ensembl_gene_id", target: "entrez_id"})
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    x.to_csv(out_path, sep="\t", index=False)
    return out_path


# id-shape signatures used to guess a matrix's namespace
_SIGNATURES = [
    ("ensembl_gene_id", lambda s: s.str.startswith("ENSCGRG").mean() > 0.5),
    ("ensembl_protein", lambda s: s.str.startswith("ENSCGRP").mean() > 0.5),
    ("refseq_mrna", lambda s: s.str.match(r"^[XN]M_\d+").mean() > 0.5),
    ("refseq_protein", lambda s: s.str.match(r"^[XN]P_\d+").mean() > 0.5),
    ("uniprot", lambda s: s.str.match(r"^[A-NR-Z0-9]{6,10}$").mean() > 0.5),
    ("entrez_id", lambda s: s.str.match(r"^\d+$").mean() > 0.5),
    ("symbol", lambda s: s.str.match(r"^(LOC\d+|[A-Za-z][\w.-]+)$").mean() > 0.5),
]


def detect_namespace(ids) -> str | None:
    """Best-guess which namespace a set of gene ids is keyed by (or None)."""
    s = pd.Series([str(i) for i in ids], dtype="string").dropna()
    if s.empty:
        return None
    for name, test in _SIGNATURES:
        try:
            if test(s):
                return name
        except Exception:
            continue
    return None


def summary() -> pd.Series:
    """Per-column non-null counts across the master table (a quick coverage view)."""
    m = load_master()
    counts = {}
    for c in NAMESPACES:
        if c in MULTI_COLS:
            counts[c] = int(m[c].notna().sum())
        else:
            counts[c] = int(m[c].notna().sum())
    counts["rows"] = len(m)
    counts["in_icho3k"] = int(m["in_icho3k"].sum())
    return pd.Series(counts)
