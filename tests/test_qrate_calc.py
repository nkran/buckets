import numpy as np
import pandas as pd

from gems import config as C
from gems import qrate_calc


def _metadata_single_vessel(gluc, lac, vcd, days=(3, 5)):
    rows = []
    for d, g, l, v in zip(days, gluc, lac, vcd):
        rows.append({
            "project": "P1", "clone": 1, "condition": "FB", "timepoint": d,
            "VCD": v, "Viability": 95, "Titer": 100.0 * d,
            "Gluc": g, "Gln": 5.0, "Glu": 1.0, "Lac": l, "NH3": 1.0,
        })
    return pd.DataFrame(rows).set_index(
        pd.Index([f"s{i}" for i in range(len(rows))], name="sample_id"))


def test_qrate_sign_and_magnitude():
    # glucose consumed 40->30 (ΔS=-10), lactate produced 2->8 (ΔS=+6), over 2 days
    meta = _metadata_single_vessel(gluc=[40, 30], lac=[2, 8], vcd=[5, 5], days=[3, 5])
    q = qrate_calc.compute_qrates(meta)
    ivcd = 0.5 * (5 + 5) * 2  # = 10
    scale = qrate_calc.FLUX_SCALE
    row = q[q["metabolite"] == "Gluc"].iloc[0]
    assert np.isclose(row["ivcd"], ivcd)
    assert np.isclose(row["q"], -10.0 / ivcd * scale)   # consumed -> negative
    lac = q[q["metabolite"] == "Lac"].iloc[0]
    assert np.isclose(lac["q"], 6.0 / ivcd * scale)      # produced -> positive
    assert row["q"] < 0 and lac["q"] > 0


def test_qrate_bounds_bracket_estimate():
    meta = _metadata_single_vessel(gluc=[40, 30], lac=[2, 8], vcd=[5, 5], days=[3, 5])
    q = qrate_calc.compute_qrates(meta)
    bounds = qrate_calc.qrate_bounds_for_sample(q, "P1", 1, "FB", 5)
    lb, ub = bounds[C.METABOLITE_EXCHANGE["Gluc"]]
    q_est = -10.0 / 10.0 * qrate_calc.FLUX_SCALE
    assert lb <= q_est <= ub
    assert lb < 0 and ub < 0  # uptake stays a net uptake band
