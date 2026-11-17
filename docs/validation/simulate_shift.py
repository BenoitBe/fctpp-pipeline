"""Simulation: calibration and power of the shift test vs the legacy rule.

    DEIMOS_DIR=/path/to/deimos python docs/validation/simulate_shift.py [--reps 3] [--out results.csv]

Model (lysate-treated TPP, per protein):
  soluble fraction f(T) = 1 / (1 + exp((T - Tm) / s)),  Tm ~ N(52, 6) °C, s ~ U(1.5, 3)
  log2 intensity = base + log2(f(T) + floor) + lysate effect + noise
  - 12 temperatures 37..67 °C, 3 conditions (Ctl, C1, C2), R lysates split in every
    condition and temperature (paired design, as in the real data)
  - lysate effect ~ N(0, 0.25) shared by the 3 conditions, drifting slowly with T
  - noise ~ N(0, 0.26) (median residual SD of the real data)
  - C1 stabilises 4 % of the proteins by ΔTm ∈ {2, 4, 6} °C; C2 has a constant +0.8
    offset on 2 % of the proteins (handling artefact: no thermal shift)
Truth for C1_vs_Ctl = the stabilised proteins. Each method is scored by sensitivity
(per ΔTm) and by its false discovery proportion; C2_vs_Ctl offsets measure how
often an artefact is called a shift.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from olympus.engine import fit_series, load_limma          # noqa: E402
from olympus.matrices import TemperatureMatrix             # noqa: E402
from olympus.shift import ShiftParams, adjacent_correlation, fits_to_long, shift_test  # noqa: E402
from olympus.temperatures import TemperatureSeries         # noqa: E402
from olympus.tpp import CallParams, FilterParams, run_tpp  # noqa: E402

TEMPS_C = np.linspace(37, 67, 12)
TEMPS = [f"T{i}" for i in range(1, 13)]
CONDS = ["C2", "C1", "Ctl"]
CONTRASTS = ["C2_vs_C1", "C2_vs_Ctl", "C1_vs_Ctl"]


def simulate(n_prot=2000, n_lysates=3, seed=0, frac_shift=0.04, frac_offset=0.02,
             noise=0.26, lysate_sd=0.25, realistic=True):
    """realistic=True adds what the Gaussian toy lacks, calibrated on the real data:
    protein-specific variances (scaled inverse chi-square, d0 = 2, as estimated by
    limma on the TPP runs), t-distributed noise (5 df), and left-censored missing
    values (below a per-temperature detection limit) imputed from a narrow normal
    under that limit (QRILC-like)."""
    rng = np.random.default_rng(seed)
    names = [f"P{i:05d}" for i in range(n_prot)]
    base = rng.normal(20, 2, n_prot)
    tm = rng.normal(52, 6, n_prot)
    slope = rng.uniform(1.5, 3.0, n_prot)
    if realistic:
        d0 = 2.0
        sd_p = noise * np.sqrt(d0 / rng.chisquare(d0, n_prot)) * 0.6   # median ≈ noise
        sd_p = np.clip(sd_p, 0.05, 2.0)
    else:
        sd_p = np.full(n_prot, noise)
    dtm = np.zeros(n_prot)
    shifted = rng.choice(n_prot, int(frac_shift * n_prot), replace=False)
    dtm[shifted] = rng.choice([2.0, 4.0, 6.0], len(shifted))
    offset = np.zeros(n_prot)
    rest = np.setdiff1d(np.arange(n_prot), shifted)
    off_idx = rng.choice(rest, int(frac_offset * n_prot), replace=False)
    offset[off_idx] = 0.8
    lys0 = rng.normal(0, lysate_sd, (n_prot, n_lysates))
    mats = {}
    for k, (t, tc) in enumerate(zip(TEMPS, TEMPS_C)):
        cols, rows = {}, []
        for L in range(n_lysates):
            lys = lys0[:, L] * (1 - 0.04 * k) + rng.normal(0, 0.05, n_prot)
            for cond in CONDS:
                tmc = tm + (dtm if cond == "C1" else 0.0)
                f = 1.0 / (1.0 + np.exp((tc - tmc) / slope))
                eps = (rng.standard_t(5, n_prot) / np.sqrt(5 / 3)) if realistic else rng.normal(0, 1, n_prot)
                y = base + np.log2(f + 0.05) + lys + sd_p * eps
                if cond == "C2":
                    y = y + offset
                lab = f"S{k + 1}_{cond}_{L + 1}"
                cols[lab] = y
                rows.append({"label": lab, "condition": cond, "replicate": str(L + 1)})
        imp = pd.DataFrame(cols, index=pd.Index(names, name="protein"))
        if realistic:
            lod = np.quantile(imp.to_numpy(), 0.03)               # ~3 % censored per temperature
            low = imp < lod
            imp = imp.mask(low, rng.normal(lod - 0.8, 0.3, imp.shape))
        mats[t] = TemperatureMatrix(t, imp, None, pd.DataFrame(rows), Path("."))
    truth = pd.DataFrame({"name": names, "dTm": dtm, "offset": offset})
    return mats, truth


def score(called: set, truth: pd.DataFrame, contrast: str) -> dict:
    if contrast == "C1_vs_Ctl":
        true = set(truth.loc[truth.dTm > 0, "name"])
        out = {"called": len(called), "false": len(called - true),
               "fdp": round(len(called - true) / max(len(called), 1), 3)}
        for d in (2.0, 4.0, 6.0):
            tt = set(truth.loc[truth.dTm == d, "name"])
            out[f"sens_dTm{int(d)}"] = round(len(called & tt) / max(len(tt), 1), 3)
        return out
    if contrast == "C2_vs_Ctl":
        off = set(truth.loc[truth.offset > 0, "name"])
        return {"called": len(called), "offset_called": len(called & off),
                "other_false": len(called - off)}
    return {"called": len(called)}


def run_once(limma, seed, n_lysates, realistic=True):
    mats, truth = simulate(seed=seed, n_lysates=n_lysates, realistic=realistic)
    temps = TemperatureSeries.from_ids(TEMPS, dict(zip(TEMPS, TEMPS_C)))
    rows = []
    fits_u = fit_series(limma, mats, TEMPS, CONTRASTS, paired=False)
    fits_p = fit_series(limma, mats, TEMPS, CONTRASTS, paired="auto")
    for label, fits in (("legacy rule, unpaired fit", fits_u), ("legacy rule, paired fit", fits_p)):
        long = fits_to_long(fits, TEMPS, CONTRASTS)
        r = run_tpp(long, temps, FilterParams(6, 6), CallParams(rule="strict"), contrasts=CONTRASTS)
        for c in ("C1_vs_Ctl", "C2_vs_Ctl"):
            rows.append({"method": label, "contrast": c, "seed": seed, "lysates": n_lysates,
                         **score(set(r.hits.loc[r.hits.contrast == c, "name"]), truth, c)})
    rho = adjacent_correlation(fits_p, TEMPS)
    for baseline in ("median", "lowest", "none"):
        for alpha in (0.05, 0.10):
            prm = ShiftParams(baseline=baseline, window=2, alpha=alpha, fc_cutoff=0.6)
            tab, _ = shift_test(fits_p, TEMPS, CONTRASTS, prm, rho)
            for c in ("C1_vs_Ctl", "C2_vs_Ctl"):
                rows.append({"method": f"shift test, baseline {baseline}, q≤{alpha}", "contrast": c,
                             "seed": seed, "lysates": n_lysates,
                             **score(set(tab.loc[(tab.contrast == c) & tab.hit, "name"]), truth, c)})
    # calibration under the null: C2_vs_C1 has no thermal signal except C1 shifts... use a pure-null
    # contrast instead: p-value uniformity on proteins with no shift and no offset (C1_vs_Ctl)
    tab, _ = shift_test(fits_p, TEMPS, ["C1_vs_Ctl"], ShiftParams(baseline="median", window=2), rho)
    nulls = set(truth.loc[(truth.dTm == 0), "name"])
    pv = tab.loc[tab.name.isin(nulls), "p_protein"].dropna()
    rows.append({"method": "calibration (null proteins, median)", "contrast": "C1_vs_Ctl", "seed": seed,
                 "lysates": n_lysates, "p<0.05": round(float((pv < 0.05).mean()), 4),
                 "p<0.01": round(float((pv < 0.01).mean()), 4), "p<0.001": round(float((pv < 0.001).mean()), 5)})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--lysates", type=int, nargs="+", default=[3, 4])
    ap.add_argument("--out", default=None)
    ap.add_argument("--deimos-dir", default=None)
    ap.add_argument("--gaussian", action="store_true", help="toy model (no heterogeneity, no censoring)")
    a = ap.parse_args()
    limma = load_limma(a.deimos_dir)
    import logging
    logging.getLogger("olympus.tpp").setLevel(logging.WARNING)
    rows = []
    for n_l in a.lysates:
        for s in range(a.reps):
            rows += run_once(limma, s, n_l, realistic=not a.gaussian)
    df = pd.DataFrame(rows)
    num = [c for c in df.columns if c not in ("method", "contrast", "seed", "lysates")]
    agg = df.groupby(["lysates", "contrast", "method"], sort=False)[num].mean().round(3)
    pd.set_option("display.width", 220)
    print(agg.dropna(axis=1, how="all").to_string())
    if a.out:
        df.to_csv(a.out, index=False)


if __name__ == "__main__":
    main()
