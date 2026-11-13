import ast

import numpy as np
import pandas as pd

from gems import gene_mapping


def test_eval_gpr_and_is_min_or_is_max():
    expr = {"g1": 10.0, "g2": 4.0, "g3": 7.0}
    # (g1 and g2) or g3  -> max(min(10,4), 7) = max(4,7) = 7
    tree = ast.parse("(g1 and g2) or g3", mode="eval")
    assert gene_mapping._eval_gpr(tree.body, expr) == 7.0
    # g1 and g2 and g3 -> min = 4
    tree2 = ast.parse("g1 and g2 and g3", mode="eval")
    assert gene_mapping._eval_gpr(tree2.body, expr) == 4.0


def test_eval_gpr_missing_gene_is_neutral():
    expr = {"g1": 5.0}
    tree = ast.parse("g1 and gX", mode="eval")   # gX absent -> ignored
    assert gene_mapping._eval_gpr(tree.body, expr) == 5.0
    tree2 = ast.parse("gX and gY", mode="eval")   # all absent -> NaN
    assert np.isnan(gene_mapping._eval_gpr(tree2.body, expr))


def test_counts_to_entrez_collapses_many_to_one():
    counts = pd.DataFrame({"s1": [1.0, 2.0, 4.0]},
                          index=["ENS1", "ENS2", "ENS3"])
    xwalk = pd.DataFrame({"ensembl_gene_id": ["ENS1", "ENS2", "ENS3"],
                          "entrez_id": ["100", "100", "200"]})
    out = gene_mapping.counts_to_entrez(counts, xwalk)
    assert out.loc["100", "s1"] == 3.0   # 1 + 2 summed
    assert out.loc["200", "s1"] == 4.0


def test_counts_to_entrez_one_ensembl_to_many_entrez():
    # real BioMart crosswalks map some Ensembl genes to several Entrez ids
    counts = pd.DataFrame({"s1": [5.0, 3.0]}, index=["ENS1", "ENS2"])
    xwalk = pd.DataFrame({"ensembl_gene_id": ["ENS1", "ENS1", "ENS2"],
                          "entrez_id": ["100", "200", "300"]})
    out = gene_mapping.counts_to_entrez(counts, xwalk)
    assert out.loc["100", "s1"] == 5.0   # ENS1 expression attributed to both
    assert out.loc["200", "s1"] == 5.0
    assert out.loc["300", "s1"] == 3.0


def test_resolve_aliases_picks_variant_carrying_flux(tmp_path):
    import json
    from gems import escher_view as ev
    mp = tmp_path / "m.json"
    mp.write_text("[]")
    (tmp_path / "m.aliases.json").write_text(json.dumps({"CS": ["CS", "CSm"]}))
    out = ev.resolve_aliases(mp, {"CS": 0.0, "CSm": -2.5, "PFK": 1.0})
    assert out["CS"] == -2.5 and out["PFK"] == 1.0
    assert ev.resolve_aliases(tmp_path / "none.json", {"CS": 0.0}) == {"CS": 0.0}
