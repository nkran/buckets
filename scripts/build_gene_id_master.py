"""Build a comprehensive gene-identifier bridge for Chinese hamster (CriGri-PICRH-1.0).

Connects, in one table, the identifier namespaces that show up across CHO transcriptomics
and metabolic modelling:

    Ensembl PICRH gene / transcript / protein  <->  NCBI Entrez GeneID  <->  gene symbol
        <->  RefSeq mRNA (NM_/XM_) & protein (NP_/XP_)  <->  UniProt  <->  iCHO3K model gene

Anchor = the **Ensembl PICRH** gene set (what real DESeq2 counts are keyed by), enriched
from Ensembl BioMart (dataset `cgpicr_gene_ensembl`) and flagged for presence in the
**iCHO3K** genome-scale model (whose genes are NCBI Entrez ids). iCHO3K genes that have no
Ensembl PICRH mapping are appended so the table is a superset of both namespaces.

Output: `data/maps/gene_id_master.tsv` (+ `.parquet`) and a short coverage summary.

Run:  python scripts/build_gene_id_master.py
"""

from __future__ import annotations

import io
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

from gems import config as C
from gems import model_setup

MARTSERVICE = "https://www.ensembl.org/biomart/martservice"
DATASET = "cgpicr_gene_ensembl"
ICHO3K_SYMBOL_MAP = C.DATA / "maps" / "icho3k_gene_id_mapping.csv"   # gene(symbol/LOC) -> entrez
OUT_TSV = C.DATA / "maps" / "gene_id_master.tsv"


CACHE_DIR = C.DATA / "maps" / "biomart_cache"


def biomart(attributes: list[str], retries: int = 3) -> pd.DataFrame:
    """Run a BioMart TSV query (cached to disk, with retries) and return a DataFrame."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache = CACHE_DIR / ("__".join(attributes)[:120] + ".tsv")
    if cache.exists() and cache.stat().st_size > 0:
        return pd.read_csv(cache, sep="\t", names=attributes, dtype=str).replace("", pd.NA)

    attrs = "".join(f'<Attribute name="{a}"/>' for a in attributes)
    query = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
             '<Query virtualSchemaName="default" formatter="TSV" header="0" '
             'uniqueRows="1" count="" datasetConfigVersion="0.6">'
             f'<Dataset name="{DATASET}" interface="default">{attrs}</Dataset></Query>')
    data = urllib.parse.urlencode({"query": query}).encode()
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(MARTSERVICE, data=data, timeout=600) as fh:
                text = fh.read().decode()
            if "ERROR" in text[:200] or "Query ERROR" in text:
                raise RuntimeError(f"BioMart error: {text[:300]}")
            cache.write_text(text)
            return pd.read_csv(io.StringIO(text), sep="\t", names=attributes, dtype=str).replace("", pd.NA)
        except Exception as exc:  # noqa: BLE001 — retry transient timeouts
            last = exc
            print(f"    biomart attempt {attempt+1}/{retries} failed: {str(exc)[:60]}")
    raise RuntimeError(f"BioMart query failed after {retries} attempts: {last}")


def _join_unique(series: pd.Series) -> str:
    vals = [v for v in dict.fromkeys(series.dropna()) if v]   # order-preserving unique
    return ";".join(vals)


def build() -> pd.DataFrame:
    print(f"Querying Ensembl BioMart ({DATASET}) …")

    # 1. gene-level attributes (one row per gene)
    gene = biomart(["ensembl_gene_id", "entrezgene_id", "entrezgene_accession",
                    "external_gene_name", "gene_biotype", "chromosome_name", "description"])
    gene = gene.rename(columns={"entrezgene_id": "entrez_id",
                                "external_gene_name": "symbol",
                                "chromosome_name": "chromosome"})
    # symbol: prefer Ensembl gene name, fall back to the NCBI accession/symbol
    gene["symbol"] = gene["symbol"].fillna(gene["entrezgene_accession"])
    gene = gene.drop(columns=["entrezgene_accession"])
    print(f"  gene-level rows: {len(gene)}  (genes: {gene.ensembl_gene_id.nunique()})")

    # 2. RefSeq (transcript-level) -> aggregate to gene.
    #    BioMart allows <=3 external-reference attributes per query, so mRNA and
    #    peptide are fetched separately (each: curated NM_/NP_ + predicted XM_/XP_).
    rs_m = biomart(["ensembl_gene_id", "refseq_mrna", "refseq_mrna_predicted"])
    rs_m["v"] = rs_m[["refseq_mrna", "refseq_mrna_predicted"]].apply(
        lambda r: ";".join(v for v in r.dropna()), axis=1)
    refseq_mrna = (rs_m.groupby("ensembl_gene_id").agg(refseq_mrna=("v", _join_unique))
                       .reset_index())
    rs_p = biomart(["ensembl_gene_id", "refseq_peptide", "refseq_peptide_predicted"])
    rs_p["v"] = rs_p[["refseq_peptide", "refseq_peptide_predicted"]].apply(
        lambda r: ";".join(v for v in r.dropna()), axis=1)
    refseq_pep = (rs_p.groupby("ensembl_gene_id").agg(refseq_protein=("v", _join_unique))
                      .reset_index())
    refseq = refseq_mrna.merge(refseq_pep, on="ensembl_gene_id", how="outer")
    print(f"  genes with any RefSeq mRNA: {(refseq.refseq_mrna.fillna('') != '').sum()}")

    # 3. Ensembl protein + UniProt (transcript/protein-level) -> aggregate to gene
    prot = biomart(["ensembl_gene_id", "ensembl_peptide_id",
                    "uniprotswissprot", "uniprotsptrembl"])
    prot["uniprot_all"] = prot[["uniprotswissprot", "uniprotsptrembl"]].apply(
        lambda r: ";".join(v for v in r.dropna()), axis=1)
    protein = (prot.groupby("ensembl_gene_id")
                   .agg(ensembl_protein=("ensembl_peptide_id", _join_unique),
                        uniprot=("uniprot_all", _join_unique))
                   .reset_index())

    master = gene.merge(refseq, on="ensembl_gene_id", how="left") \
                 .merge(protein, on="ensembl_gene_id", how="left")

    # 4. iCHO3K model: mark genes present + reaction counts + the model's own symbol map
    model = model_setup.get_base_model()
    model_entrez = {g.id for g in model.genes}
    rxn_count = {g.id: len(g.reactions) for g in model.genes}
    master["entrez_id"] = master["entrez_id"].astype("string")
    master["in_icho3k"] = master["entrez_id"].isin(model_entrez)
    master["icho3k_n_reactions"] = master["entrez_id"].map(rxn_count).astype("Int64")

    # 5. append iCHO3K model genes that have NO Ensembl PICRH mapping (so the table is a
    #    superset of both namespaces), backfilling symbol from the model's gene_id_mapping.
    sym_map = pd.read_csv(ICHO3K_SYMBOL_MAP, dtype=str)   # columns: gene, entrez_gene_id
    sym_map["entrez_gene_id"] = sym_map["entrez_gene_id"].str.split(".").str[0]
    entrez2sym = dict(zip(sym_map["entrez_gene_id"], sym_map["gene"]))
    mapped_entrez = set(master.loc[master["in_icho3k"], "entrez_id"].dropna())
    missing = sorted(model_entrez - mapped_entrez)
    if missing:
        extra = pd.DataFrame({
            "ensembl_gene_id": pd.NA, "entrez_id": missing,
            "symbol": [entrez2sym.get(e, pd.NA) for e in missing],
            "gene_biotype": "protein_coding", "in_icho3k": True,
            "icho3k_n_reactions": [rxn_count.get(e) for e in missing],
        })
        master = pd.concat([master, extra], ignore_index=True)
    print(f"  iCHO3K model genes: {len(model_entrez)}  "
          f"(reachable via Ensembl PICRH: {len(mapped_entrez)}, "
          f"appended without Ensembl: {len(missing)})")

    cols = ["ensembl_gene_id", "entrez_id", "symbol", "gene_biotype", "chromosome",
            "refseq_mrna", "refseq_protein", "ensembl_protein", "uniprot",
            "in_icho3k", "icho3k_n_reactions", "description"]
    return master[[c for c in cols if c in master.columns]]


def main() -> None:
    master = build()
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    master.to_csv(OUT_TSV, sep="\t", index=False)
    try:
        master.to_parquet(OUT_TSV.with_suffix(".parquet"), index=False)
    except Exception:
        pass

    n = len(master)
    print(f"\nWrote {OUT_TSV}  ({n} rows)")
    print("── coverage summary ──")
    print(f"  rows (genes)                 : {n}")
    print(f"  with Ensembl PICRH id        : {master.ensembl_gene_id.notna().sum()}")
    print(f"  with Entrez id               : {master.entrez_id.notna().sum()}")
    print(f"  with a gene symbol           : {master.symbol.notna().sum()}")
    print(f"  with RefSeq mRNA             : {(master.get('refseq_mrna').fillna('') != '').sum()}")
    print(f"  with UniProt                 : {(master.get('uniprot').fillna('') != '').sum()}")
    print(f"  in iCHO3K model              : {int(master.in_icho3k.sum())}")
    inm = master[master.in_icho3k]
    print(f"    …of which with Ensembl id  : {inm.ensembl_gene_id.notna().sum()} "
          f"({inm.ensembl_gene_id.notna().mean():.1%})")
    print(f"    …of which with RefSeq mRNA : {(inm.get('refseq_mrna').fillna('') != '').sum()}")


if __name__ == "__main__":
    main()
