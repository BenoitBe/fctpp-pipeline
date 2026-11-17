import numpy as np
import pandas as pd
import pytest

from olympus.engine import complete_blocks, fit_series, fit_temperature, parse_contrast
from olympus.shift import (ShiftParams, adjacent_correlation, bh, fits_to_long, permutation_null,
                           shift_test, simes, swap_patterns)
from olympus.synthetic import simulate_intensity_series

TEMPS = [f"T{i}" for i in range(1, 13)]
CONTRASTS = ["C2_vs_C1", "C2_vs_Ctl", "C1_vs_Ctl"]


# ---------------------------------------------------------------- pure helpers

def test_bh_matches_p_adjust():
    p = np.array([0.01, 0.04, 0.03, 0.2, np.nan])
    q = bh(p)                                   # R: p.adjust(c(.01,.04,.03,.2), "BH")
    np.testing.assert_allclose(q[:4], [0.04, 0.05333333, 0.05333333, 0.2], rtol=1e-6)
    assert np.isnan(q[4])


def test_simes():
    s = simes(np.array([[0.01, 0.5, 0.9], [np.nan, np.nan, np.nan], [0.02, 0.03, np.nan]]))
    assert s[0] == pytest.approx(0.03)          # min(0.01*3/1, 0.5*3/2, 0.9)
    assert np.isnan(s[1])
    assert s[2] == pytest.approx(0.03)          # min(0.02*2, 0.03*2/2)


def test_swap_patterns_dedupe_complements():
    assert len(swap_patterns(["1", "2", "3"])) == 3      # 6 non-trivial subsets / complement pairs
    assert len(swap_patterns(["1", "2", "3", "4"])) == 7
    assert len(swap_patterns([str(i) for i in range(8)], max_patterns=10)) == 10


def test_complete_blocks_and_contrast_parsing():
    assert complete_blocks(["A", "B", "A", "B"], ["1", "1", "2", "2"])
    assert not complete_blocks(["A", "B", "A", "B"], ["1", "2", "3", "4"])
    assert parse_contrast("C1_vs_Ctl") == ("C1", "Ctl")
    with pytest.raises(ValueError):
        parse_contrast("C1-Ctl")


def test_params_validation():
    with pytest.raises(ValueError):
        ShiftParams(baseline="mean").validate()
    with pytest.raises(ValueError):
        ShiftParams(localisation_p=0).validate()


# ---------------------------------------------------------------- with Deimos' engine

def test_paired_fit_auto(limma_engine, sim_series):
    mats, _ = sim_series
    tm = mats["T5"]
    f = fit_temperature(limma_engine, "T5", tm.imputed, tm.samples["condition"].tolist(),
                        tm.block_labels(), CONTRASTS, paired="auto")
    u = fit_temperature(limma_engine, "T5", tm.imputed, tm.samples["condition"].tolist(),
                        tm.block_labels(), CONTRASTS, paired=False)
    assert f.paired and not u.paired
    np.testing.assert_allclose(f.diff, u.diff, atol=1e-10)     # balanced blocks: same estimates
    assert np.median(f.se) < np.median(u.se)                   # lysate effect removed


def test_shift_test_detects_shifts_and_ignores_offsets(limma_engine, sim_series):
    mats, truth = sim_series
    fits = fit_series(limma_engine, mats, TEMPS, CONTRASTS, paired="auto")
    rho = adjacent_correlation(fits, TEMPS)
    assert all(0.0 <= v <= 0.95 for v in rho.values())
    tab, wins = shift_test(fits, TEMPS, CONTRASTS, ShiftParams(), rho)
    hits = set(tab.loc[(tab.contrast == "C1_vs_Ctl") & tab.hit, "name"])
    big = set(truth.loc[truth.dTm >= 4, "name"])
    assert len(hits & big) / len(big) >= 0.7
    fdp = len(hits - set(truth.loc[truth.dTm > 0, "name"])) / max(len(hits), 1)
    assert fdp <= 0.15
    off = set(truth.loc[truth.offset > 0, "name"])
    c2 = set(tab.loc[(tab.contrast == "C2_vs_Ctl") & tab.hit, "name"])
    assert len(c2 & off) <= 0.25 * len(off)                     # offsets are (mostly) not shifts
    # without the localisation check, offsets come back
    tab2, _ = shift_test(fits, TEMPS, ["C2_vs_Ctl"], ShiftParams(localisation_p=None), rho)
    assert len(set(tab2.loc[tab2.hit, "name"]) & off) >= len(off) // 2
    assert {"window", "dd", "se", "df", "p"} <= set(wins.columns)


def test_shift_test_null_calibration(limma_engine):
    mats, _ = simulate_intensity_series(n_prot=800, seed=5, frac_shift=0.0, frac_offset=0.0)
    fits = fit_series(limma_engine, mats, TEMPS, CONTRASTS, paired="auto")
    tab, _ = shift_test(fits, TEMPS, CONTRASTS, ShiftParams(), adjacent_correlation(fits, TEMPS))
    assert int(tab.hit.sum()) <= 2
    assert (tab.p_protein.dropna() < 0.05).mean() < 0.08


def test_permutation_null_counts(limma_engine, sim_series):
    mats, _ = sim_series
    calls = []

    def evaluate(f, c):
        calls.append(c)
        long = fits_to_long(f, TEMPS, [c])
        assert {f"{c}_diff", f"{c}_p.val", f"{c}_p.adj"} <= set(long.columns)
        return 1

    out = permutation_null(limma_engine, mats, TEMPS, ["C1_vs_Ctl"], evaluate)
    assert out == {"C1_vs_Ctl": [1, 1, 1]}                      # 3 lysates -> 3 patterns


# ---------------------------------------------------------------- pairing, windows, shared component

def test_block_correlation_decides_pairing(sim_series):
    from olympus.engine import block_correlation
    mats, _ = sim_series
    bc = block_correlation(mats, TEMPS)
    assert bc.complete and bc.n_blocks == 3
    assert bc.rho > 0.2 and bc.p < 1e-4 and bc.decide()
    assert bc.different < 0                                     # ≈ −ρ/(n−1)
    from dataclasses import replace
    rng = np.random.default_rng(3)                              # pure noise: no lysate effect
    m0 = {t: replace(tm, imputed=pd.DataFrame(rng.normal(20, 0.3, tm.imputed.shape), index=tm.imputed.index,
                                              columns=tm.imputed.columns)) for t, tm in mats.items()}
    b0 = block_correlation(m0, TEMPS)
    assert abs(b0.rho) < 0.03 and abs(b0.delta.mean()) < 0.03


def test_window_sizes_and_override():
    assert ShiftParams().window_sizes == (2, 3, 4)
    assert ShiftParams(window=3).window_sizes == (3,)           # single size overrides windows
    assert ShiftParams(windows=(4, 2)).window_sizes == (2, 4)
    with pytest.raises(ValueError):
        ShiftParams(windows=(0, 2)).validate()


def test_multi_window_scan(limma_engine, sim_series):
    mats, truth = sim_series
    fits = fit_series(limma_engine, mats, TEMPS, ["C1_vs_Ctl"], paired=True)
    rho = adjacent_correlation(fits, TEMPS)
    t1, w1 = shift_test(fits, TEMPS, ["C1_vs_Ctl"], ShiftParams(window=2), rho)
    t3, w3 = shift_test(fits, TEMPS, ["C1_vs_Ctl"], ShiftParams(), rho)
    assert set(w3["size"]) == {2, 3, 4} and set(w1["size"]) == {2}
    assert (t1["dd_peak"] == t1["dd_best"]).all()               # single size: peak = best window
    sub = t3.dropna(subset=["dd_best"])
    assert (sub["dd_peak"].abs() >= sub["dd_best"].abs() - 1e-12).all()
    true = set(truth.loc[truth.dTm > 0, "name"])
    h1, h3 = set(t1.loc[t1.hit, "name"]), set(t3.loc[t3.hit, "name"])
    assert len(h3 & true) >= len(h1 & true) - 2                 # no loss of power
    assert len(h3 - true) / max(len(h3), 1) <= 0.15


def test_shared_component_flags_common_pattern(limma_engine):
    from olympus.shared import SharedParams, delta_matrix, shared_component
    mats, truth = simulate_intensity_series(n_prot=800, seed=7, shared_frac=0.08)
    fits = fit_series(limma_engine, mats, TEMPS, ["C1_vs_Ctl"], paired=True)
    names = sorted(mats["T1"].imputed.index)
    sc = shared_component(delta_matrix(fits, TEMPS, "C1_vs_Ctl", names), "C1_vs_Ctl", SharedParams(n_null=50))
    assert sc.strong and sc.variance_explained > 1.5 * sc.null_q95
    assert sc.loadings.idxmax() in ("T5", "T6", "T11")
    flags = sc.flags(0.7)
    shared = set(truth.loc[truth.shared, "name"])
    shifted = set(truth.loc[truth.dTm > 0, "name"])
    assert flags.reindex(list(shared)).mean() >= 0.6
    assert flags.reindex(list(shifted)).mean() <= 0.25
    # without the planted pattern: no strong component, nothing flagged
    m0, _ = simulate_intensity_series(n_prot=800, seed=7)
    f0 = fit_series(limma_engine, m0, TEMPS, ["C1_vs_Ctl"], paired=True)
    s0 = shared_component(delta_matrix(f0, TEMPS, "C1_vs_Ctl", names), "C1_vs_Ctl", SharedParams(n_null=50))
    assert not s0.strong and not s0.flags(0.7).any()
