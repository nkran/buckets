"""End-to-end smoke test on synthetic data (small sample subset for speed)."""

import pandas as pd
import pytest

from gems import config as C
from gems import pipeline


@pytest.mark.slow
def test_pipeline_end_to_end_small(tmp_path):
    # everything goes to a temp dir: tests must not touch data/ or results/
    out = pipeline.run(synthetic=True, limit=6, make_figures=True,
                       out_dir=tmp_path / "results", synth_dir=tmp_path / "synthetic")
    assert (tmp_path / "results" / "fluxes" / "flux_table.csv").exists()
    ft = out["flux_table"]
    assert isinstance(ft, pd.DataFrame)
    assert not ft.empty
    # biomass and at least one TCA reaction present in the panel
    assert (ft["reaction_id"] == C.BIOMASS_RXN).any()
    assert (ft["subsystem"] == "CITRIC ACID CYCLE").any()
    # most samples should be feasible
    assert out["manifest"]["n_infeasible"] <= out["manifest"]["n_samples"] // 2
    # figures written
    assert len(out["figures"]) >= 2
