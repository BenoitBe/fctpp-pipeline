import numpy as np
import pandas as pd
import pytest

from fctpp.inputs import discover_deimos_runs, load_deimos_runs
from fctpp.temperatures import TemperatureSeries
from fctpp.tpp import (CallParams, FilterParams, call_protein, detect_contrasts, longest_run,
                       presence_stats, run_tpp)


# --------------------------------------------------------------------------- helpers

def r_tidy_tpp_threshold(pos, diff, logp, fc, lp):
    """Literal transcription of the R loop (tidy_tpp.R, lines 989-1000)."""
    cond_k = (np.abs(diff) >= fc) & (logp >= lp)
    cond_suiv = np.append(cond_k[1:], False)
    return bool(np.any(cond_k & cond_suiv))


def test_longest_run_matches_r_rle():
    assert longest_run(np.array([], int)) == 0
    assert longest_run(np.array([3])) == 1                    # R: 1 when no diff == 1
    assert longest_run(np.array([0, 1, 2, 4, 5])) == 3
    assert longest_run(np.array([0, 2, 4])) == 1
    assert longest_run(np.arange(12)) == 12


def test_detect_contrasts_ignores_secondary_columns():
    cols = ["name", "A_vs_B_diff", "A_vs_B_p.val", "A_vs_B_p.adj", "A_vs_B_p.val_limma",
            "C_vs_B_p.val", "Pi_Score_A_vs_B"]
    assert detect_contrasts(cols) == ["A_vs_B"]               # C_vs_B has no _diff


@pytest.mark.parametrize("seed", range(40))
def test_r_legacy_rule_equals_r_loop(seed):
    rng = np.random.default_rng(seed)
    n = rng.integers(2, 12)
    pos = np.sort(rng.choice(12, size=n, replace=False))
    diff = rng.normal(0, 0.8, n)
    logp = rng.exponential(1.5, n)
    p = CallParams(rule="r_legacy")
    got = call_protein(pos, diff, logp, p)["call"] == "Over"
    assert got == r_tidy_tpp_threshold(pos, diff, logp, 0.6, 2.0)


def test_strict_rule_rejects_gap_and_sign_flip():
    p = CallParams(rule="strict")
    # gap: positions 4 and 6 are successive rows but not adjacent temperatures
    r = call_protein(np.array([3, 4, 6, 7]), np.array([0.0, 1.0, 1.0, 0.0]),
                     np.array([0.0, 3.0, 3.0, 0.0]), p)
    assert r["call"] == "Under" and r["over_r_legacy"] and not r["over_strict"]
    # sign flip on adjacent temperatures
    r = call_protein(np.array([3, 4]), np.array([0.8, -0.8]), np.array([3.0, 3.0]), p)
    assert r["call"] == "Under" and r["over_r_legacy"]
    # genuine run
    r = call_protein(np.array([5, 6, 7]), np.array([-0.9, -1.1, -0.7]), np.array([2.5, 4, 2.1]), p)
    assert r["call"] == "Over" and r["direction"] == "down" and r["run_length"] == 3


def test_min_consecutive_and_nan():
    p = CallParams(rule="strict", min_consecutive=3)
    r = call_protein(np.array([0, 1, 2]), np.array([1.0, 1.0, np.nan]), np.array([3.0, 3.0, 3.0]), p)
    assert r["call"] == "Under" and r["n_sig_temps"] == 2
    r = call_protein(np.array([0, 1, 2]), np.array([1.0, 1.0, 1.0]), np.array([3.0, 3.0, 3.0]), p)
    assert r["call"] == "Over"


def test_callparams_validation():
    with pytest.raises(ValueError):
        CallParams(rule="loose").validate()
    with pytest.raises(ValueError):
        CallParams(p_column="q").validate()


# --------------------------------------------------------------------------- run_tpp

@pytest.fixture(scope="module")
def loaded(deimos_root):
    return load_deimos_runs(discover_deimos_runs(deimos_root))


def _hits(res, contrast):
    return set(res.hits.loc[res.hits["contrast"] == contrast, "name"])


def test_run_tpp_strict_recovers_planted_shifts(loaded):
    res = run_tpp(loaded.results, contrasts=loaded.contrasts)
    assert res.temperatures.ids[:3] == ("T1", "T2", "T3") and res.temperatures.ids[-1] == "T12"
    assert _hits(res, "C1_vs_Ctl") == {"P0001", "P0002", "P0003", "P0004", "P0005", "P0009", "P0013"}
    assert _hits(res, "C2_vs_Ctl") == {"P0006", "P0007", "P0008"}
    h = res.hits.set_index(["contrast", "name"])
    assert h.loc[("C1_vs_Ctl", "P0009"), "shift_class"] == "abundance-like"
    assert h.loc[("C1_vs_Ctl", "P0001"), "shift_class"] == "thermal"
    assert h.loc[("C1_vs_Ctl", "P0001"), "run_start"] == "T7"
    assert h.loc[("C1_vs_Ctl", "P0001"), "run_end"] == "T10"
    assert h.loc[("C2_vs_Ctl", "P0006"), "direction"] == "down"
    assert h.loc[("C1_vs_Ctl", "P0013"), "imputed_in_run"] is True
    assert h.loc[("C1_vs_Ctl", "P0001"), "imputed_in_run"] is False
    # P0012 is quantified at 5 temperatures only
    st = res.stats.set_index("name")
    assert not st.loc["P0012", "retained"] and st.loc["P0012", "nb_temp"] == 5
    assert "P0012" not in set(res.data["name"])
    assert res.n_lost == 1


def test_run_tpp_r_legacy_adds_gap_and_flip_artefacts(loaded):
    res = run_tpp(loaded.results, call_params=CallParams(rule="r_legacy"), contrasts=loaded.contrasts)
    extra = _hits(res, "C1_vs_Ctl") - {"P0001", "P0002", "P0003", "P0004", "P0005", "P0009", "P0013"}
    assert extra == {"P0010", "P0011"}
    s = res.summary.set_index("contrast")
    assert s.loc["C1_vs_Ctl", "hits_r_legacy_rule"] == 9
    assert s.loc["C1_vs_Ctl", "hits_strict_rule"] == 7


def test_threshold_columns_are_r_compatible(loaded):
    res = run_tpp(loaded.results, contrasts=loaded.contrasts)
    for c in loaded.contrasts:
        for s in ("diff", "p.adj", "p.val", "logpval", "threshold"):
            assert f"{c}_{s}" in res.data.columns
        assert set(res.data[f"{c}_threshold"]) <= {"Over", "Under"}
    over = res.data.loc[res.data["C1_vs_Ctl_threshold"] == "Over", "name"].unique()
    assert set(over) == _hits(res, "C1_vs_Ctl")
    np.testing.assert_allclose(res.data["C1_vs_Ctl_logpval"], -np.log10(res.data["C1_vs_Ctl_p.val"]))


def test_presence_filter_parameters(loaded):
    res = run_tpp(loaded.results, filter_params=FilterParams(min_temperatures=12, min_streak=12),
                  contrasts=loaded.contrasts)
    assert {"P0011", "P0012"}.isdisjoint(set(res.data["name"]))


def test_errors_on_bad_input():
    df = pd.DataFrame({"Temperature": ["T1", "T1"], "name": ["A", "A"],
                       "X_vs_Y_diff": [1, 1], "X_vs_Y_p.val": [0.1, 0.1]})
    with pytest.raises(ValueError, match="duplicated"):
        run_tpp(df)
    with pytest.raises(ValueError, match="not in the series"):
        run_tpp(df.iloc[:1], temps=TemperatureSeries.from_ids(["T2"]))
    with pytest.raises(ValueError, match="No contrast"):
        run_tpp(pd.DataFrame({"Temperature": ["T1"], "name": ["A"]}))


def test_presence_stats_uses_series_adjacency():
    temps = TemperatureSeries.from_ids(["T1", "T2", "T3", "T5"])   # T4 absent from the experiment
    df = pd.DataFrame({"Temperature": ["T3", "T5"], "name": ["A", "A"]})
    st = presence_stats(df, temps)
    assert st.loc[0, "max_consecutive"] == 2        # T3 and T5 are adjacent in this series
