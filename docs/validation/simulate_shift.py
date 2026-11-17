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
from olympus.engine import block_correlation, fit_series, load_limma  # noqa: E402
from olympus.synthetic import SIM_TEMPS_C, simulate_intensity_series  # noqa: E402
from olympus.shift import ShiftParams, adjacent_correlation, fits_to_long, shift_test  # noqa: E402
from olympus.temperatures import TemperatureSeries         # noqa: E402
from olympus.tpp import CallParams, FilterParams, run_tpp  # noqa: E402

TEMPS = [f"T{i}" for i in range(1, 13)]
TEMPS_C = SIM_TEMPS_C
CONTRASTS = ["C2_vs_C1", "C2_vs_Ctl", "C1_vs_Ctl"]


def simulate(n_prot=2000, n_lysates=3, seed=0, realistic=True):
    return simulate_intensity_series(n_prot=n_prot, n_lysates=n_lysates, seed=seed, realistic=realistic)


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
    bc = block_correlation(mats, TEMPS)
    rows.append({"method": "pairing decision (auto)", "contrast": "-", "seed": seed, "lysates": n_lysates,
                 "paired": int(bc.decide()), "block_rho": round(bc.rho, 3), "block_p": bc.p})
    for wins in ((2, 3), (2, 3, 4)):
        prm = ShiftParams(baseline="none", windows=wins, alpha=0.05, fc_cutoff=0.6, localisation_p=0.01)
        tab, _ = shift_test(fits_p, TEMPS, CONTRASTS, prm, rho)
        for c in ("C1_vs_Ctl", "C2_vs_Ctl"):
            rows.append({"method": f"shift test, baseline none + localisation p≤0.01, windows "
                                   f"{'/'.join(map(str, wins))}, q≤0.05", "contrast": c,
                         "seed": seed, "lysates": n_lysates,
                         **score(set(tab.loc[(tab.contrast == c) & tab.hit, "name"]), truth, c)})
    for baseline, loc in (("none", 0.01), ("none", None), ("median", None), ("lowest", None)):
        for alpha in (0.05, 0.10):
            prm = ShiftParams(baseline=baseline, window=2, alpha=alpha, fc_cutoff=0.6, localisation_p=loc)
            tab, _ = shift_test(fits_p, TEMPS, CONTRASTS, prm, rho)
            for c in ("C1_vs_Ctl", "C2_vs_Ctl"):
                name = f"shift test, baseline {baseline}" + (f" + localisation p≤{loc}" if loc else "")
                rows.append({"method": f"{name}, q≤{alpha}", "contrast": c,
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
