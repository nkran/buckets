"""Map Ensembl PICRH gene expression onto iCHO3K reaction scores.

Two hops then GPR resolution:
  1. Ensembl PICRH (``ENSCGRG…``) -> NCBI Entrez, via a cached crosswalk table
     (built from Ensembl BioMart / PICRH annotation for real data; generated
     alongside the synthetic dataset otherwise).
  2. Entrez expression -> per-reaction score through the model's GPR, resolving
     AND with ``min`` and OR with ``max`` (standard GIMME-style continuous mapping).
"""

from __future__ import annotations

import ast
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C


@dataclass
class MappingReport:
    n_expressed_genes: int = 0
    n_mapped_to_entrez: int = 0
    n_model_genes: int = 0
    n_model_genes_with_data: int = 0
    n_reactions_scored: int = 0
    details: dict = field(default_factory=dict)

    @property
    def entrez_coverage(self) -> float:
        return self.n_mapped_to_entrez / max(self.n_expressed_genes, 1)

    @property
    def model_gene_coverage(self) -> float:
        return self.n_model_genes_with_data / max(self.n_model_genes, 1)

    def as_dict(self) -> dict:
        return {
            "n_expressed_genes": self.n_expressed_genes,
            "n_mapped_to_entrez": self.n_mapped_to_entrez,
            "entrez_coverage": round(self.entrez_coverage, 4),
            "n_model_genes": self.n_model_genes,
            "n_model_genes_with_data": self.n_model_genes_with_data,
            "model_gene_coverage": round(self.model_gene_coverage, 4),
            "n_reactions_scored": self.n_reactions_scored,
        }


def load_crosswalk(path: Path = C.CROSSWALK_TSV) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    if not {"ensembl_gene_id", "entrez_id"} <= set(df.columns):
        raise ValueError("crosswalk must have columns ensembl_gene_id, entrez_id")
    return df.dropna(subset=["ensembl_gene_id", "entrez_id"]).drop_duplicates()


def counts_to_entrez(counts: pd.DataFrame, crosswalk: pd.DataFrame) -> pd.DataFrame:
    """Collapse Ensembl-indexed counts to Entrez-indexed.

    Handles both directions of a real crosswalk: many Ensembl -> one Entrez are
    **summed**, and one Ensembl -> many Entrez attributes that gene's expression to
    each Entrez id.
    """
    xref = crosswalk[["ensembl_gene_id", "entrez_id"]].dropna().drop_duplicates()
    value_cols = list(counts.columns)
    merged = xref.merge(counts, left_on="ensembl_gene_id", right_index=True)
    if merged.empty:
        return pd.DataFrame(columns=value_cols)
    return merged.groupby("entrez_id")[value_cols].sum()


# ---------------------------------------------------------------------------
# GPR resolution
# ---------------------------------------------------------------------------
def _eval_gpr(node, expr: dict[str, float]) -> float:
    """Recursively evaluate a GPR AST node into a continuous score.

    AND -> min, OR -> max over gene scores. Missing genes are ignored within an
    expression; a node reducing to no data returns NaN (neutral / unmapped).
    """
    if isinstance(node, ast.Expression):
        return _eval_gpr(node.body, expr)
    if isinstance(node, ast.Name):
        return expr.get(node.id, np.nan)
    if isinstance(node, ast.BoolOp):
        vals = [_eval_gpr(v, expr) for v in node.values]
        vals = [v for v in vals if not np.isnan(v)]
        if not vals:
            return np.nan
        return min(vals) if isinstance(node.op, ast.And) else max(vals)
    # unexpected node -> neutral
    return np.nan


def reaction_scores(model, entrez_expr: pd.Series) -> dict[str, float]:
    """Score every reaction from a per-Entrez expression Series.

    Uses cobra's parsed GPR AST (``reaction.gpr.body``), whose ``ast.Name`` nodes
    carry the Entrez gene ids directly (these are digit-leading and so cannot be
    re-parsed with ``ast.parse``). Reactions without a GPR score to NaN.
    """
    expr = {str(k): float(v) for k, v in entrez_expr.items()}
    scores: dict[str, float] = {}
    for rxn in model.reactions:
        if not rxn.gene_reaction_rule:
            scores[rxn.id] = np.nan
            continue
        body = getattr(rxn.gpr, "body", None)
        if body is None:
            scores[rxn.id] = np.nan
            continue
        scores[rxn.id] = _eval_gpr(body, expr)
    return scores


def map_expression(model, counts: pd.DataFrame, crosswalk: pd.DataFrame,
                   sample_id: str) -> tuple[dict[str, float], MappingReport]:
    """Map one sample's expression to reaction scores + a coverage report."""
    entrez_expr = counts_to_entrez(counts[[sample_id]], crosswalk)[sample_id]
    scores = reaction_scores(model, entrez_expr)

    model_genes = {g.id for g in model.genes}
    # count expressed genes that mapped (not Entrez ids), so coverage stays ≤ 1
    xref = crosswalk[["ensembl_gene_id", "entrez_id"]].dropna()
    n_expressed_mapped = int(counts.index.isin(set(xref["ensembl_gene_id"])).sum())
    rep = MappingReport(
        n_expressed_genes=int(counts.shape[0]),
        n_mapped_to_entrez=n_expressed_mapped,
        n_model_genes=len(model_genes),
        n_model_genes_with_data=len(model_genes & set(entrez_expr.index.astype(str))),
        n_reactions_scored=int(sum(1 for v in scores.values() if not np.isnan(v))),
    )
    return scores, rep
