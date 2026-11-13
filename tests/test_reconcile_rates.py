"""reconcile_rates: measured exchange rates -> closest feasible bounds (toy model)."""

import cobra
import pytest

from gems import qrate_calc


def _toy_model():
    # A_e -> B_e, with biomass draining a little of A. EX_A can only be taken up.
    m = cobra.Model("toy")
    a, b = cobra.Metabolite("a_e", compartment="e"), cobra.Metabolite("b_e", compartment="e")
    conv = cobra.Reaction("CONV", lower_bound=0, upper_bound=1000)
    conv.add_metabolites({a: -1, b: 1})
    bio = cobra.Reaction("BIO", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({a: -1})
    ex_a = cobra.Reaction("EX_a_e", lower_bound=-10, upper_bound=0)
    ex_a.add_metabolites({a: -1})
    ex_b = cobra.Reaction("EX_b_e", lower_bound=0, upper_bound=1000)
    ex_b.add_metabolites({b: -1})
    m.add_reactions([conv, bio, ex_a, ex_b])
    m.objective = "BIO"
    return m


def test_consistent_rates_are_kept():
    m = _toy_model()
    bounds, rep = qrate_calc.reconcile_rates(
        m, {"EX_a_e": -2.0, "EX_b_e": 1.5}, min_growth=0.1, biomass_rxn="BIO")
    assert rep["deviation"].abs().sum() == pytest.approx(0, abs=1e-6)
    lb, ub = bounds["EX_a_e"]
    assert lb <= -2.0 <= ub


def test_inconsistent_rates_are_minimally_adjusted():
    m = _toy_model()
    # secreting 3 of B from 2 of A is impossible: total adjustment must be exactly 1
    bounds, rep = qrate_calc.reconcile_rates(
        m, {"EX_a_e": -2.0, "EX_b_e": 3.0}, min_growth=0.0, biomass_rxn="BIO")
    assert rep["deviation"].abs().sum() == pytest.approx(1.0, abs=1e-6)
    with m:
        for rid, b in bounds.items():
            m.reactions.get_by_id(rid).bounds = b
        assert m.slim_optimize() == m.slim_optimize()   # feasible (not NaN)


def test_returns_none_when_growth_is_impossible():
    m = _toy_model()
    m.reactions.EX_a_e.upper_bound = 0
    out = qrate_calc.reconcile_rates(m, {"EX_b_e": 0.0}, min_growth=50.0, biomass_rxn="BIO")
    assert out == (None, None)
