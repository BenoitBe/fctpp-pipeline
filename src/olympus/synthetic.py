"""Synthetic data with planted thermal shifts (tests, CI smoke test, demo).

``make_deimos_like_runs`` writes one Deimos-like workbook per temperature
(``Differential_Expression`` + ``Methods`` sheets only). The planted cases cover
every branch of the call logic:

=========  ==========================================================  =======================
protein    pattern (contrast C1_vs_Ctl unless stated)                  expected (strict rule)
=========  ==========================================================  =======================
P0001-05   +1.0 log2FC, p=1e-4 at T7..T10                              Over, up, thermal
P0006-08   -0.9 at T5..T7 in C2_vs_Ctl                                 Over, down, thermal
P0009      +1.2 at every temperature                                   Over, abundance-like
P0010      +0.8 at T4 then -0.8 at T5 (sign flip)                      Under (r_legacy: Over)
P0011      absent at T6, +0.9 at T5 and T7 (gap)                       Under (r_legacy: Over)
P0012      quantified at 5 temperatures only                           filtered out
P0013      as P0001 but num_NAs = 2 at T8                              Over, imputed_in_run
=========  ==========================================================  =======================
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

TEMPS = [f"T{i}" for i in range(1, 13)]
CELSIUS = dict(zip(TEMPS, [37.0, 40.4, 44.0, 46.9, 49.8, 52.9, 55.5, 58.6, 62.0, 66.3, 67.7, 70.0]))
CONTRASTS = ["C2_vs_C1", "C2_vs_Ctl", "C1_vs_Ctl"]


def _planted(n_proteins: int, seed: int) -> dict[str, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    names = [f"P{i:04d}" for i in range(1, n_proteins + 1)]
    genes = [f"gen{i:04d}" for i in range(1, n_proteins + 1)]
    tables = {}
    for ti, t in enumerate(TEMPS, start=1):
        d = {c: rng.normal(0, 0.15, n_proteins) for c in ("C1_vs_Ctl", "C2_vs_Ctl")}
        p = {c: rng.uniform(0.02, 1, n_proteins) for c in ("C1_vs_Ctl", "C2_vs_Ctl")}
        nn = np.zeros(n_proteins, int)

        def put(idx, c, fc, pv):
            d[c][idx] = fc
            p[c][idx] = pv

        if 7 <= ti <= 10:
            for i in range(0, 5):
                put(i, "C1_vs_Ctl", 1.0, 1e-4)
            put(12, "C1_vs_Ctl", 1.0, 1e-4)
            if ti == 8:
                nn[12] = 2
        if 5 <= ti <= 7:
            for i in range(5, 8):
                put(i, "C2_vs_Ctl", -0.9, 1e-4)
        put(8, "C1_vs_Ctl", 1.2, 1e-4)
        if ti == 4:
            put(9, "C1_vs_Ctl", 0.8, 1e-3)
        if ti == 5:
            put(9, "C1_vs_Ctl", -0.8, 1e-3)
        if ti in (5, 7):
            put(10, "C1_vs_Ctl", 0.9, 1e-3)

        diff_c2c1 = d["C2_vs_Ctl"] - d["C1_vs_Ctl"]
        p_c2c1 = np.minimum(1, np.maximum(rng.uniform(0.02, 1, n_proteins),
                                          np.where(np.abs(diff_c2c1) > 0.7, 1e-3, 0)))
        df = pd.DataFrame({
            "name": names, "Protein.Group": names, "Genes": genes,
            "First.Protein.Description": [f"Synthetic protein {g}" for g in genes],
            "imputed": np.where(nn > 0, "VRAI", "FAUX"), "num_NAs": nn,
        })
        for c, dd, pp in (("C2_vs_C1", diff_c2c1, p_c2c1),
                          ("C2_vs_Ctl", d["C2_vs_Ctl"], p["C2_vs_Ctl"]),
                          ("C1_vs_Ctl", d["C1_vs_Ctl"], p["C1_vs_Ctl"])):
            df[f"{c}_diff"] = dd
            df[f"{c}_p.val"] = pp
            df[f"{c}_p.adj"] = np.minimum(1, pp * 20)
            df[f"Pi_Score_{c}"] = np.abs(dd) * -np.log10(pp)
        drop = []
        if ti == 6:
            drop.append(10)                 # P0011 gap at T6
        if ti > 5:
            drop.append(11)                 # P0012 only T1..T5
        tables[t] = df.drop(index=drop).reset_index(drop=True)
    return tables


def make_deimos_like_runs(root: str | Path, n_proteins: int = 300, seed: int = 7,
                          temps: list[str] | None = None, layout: str = "nested") -> Path:
    """One Deimos-like workbook per temperature (default T1..T12).

    layout: ``nested`` → ``<root>/<T>/proteomics_output/ProteomicAnalysis_Results.xlsx``
    (Deimos default out_dir inside one folder per temperature); ``flat`` →
    ``<root>/<T>_ProteomicAnalysis_Results.xlsx``.
    """
    root = Path(root)
    methods = pd.DataFrame(
        [["LFQ PROTEOMICS PIPELINE — Methods & Parameters", None, None],
         ["2. DATA PREPROCESSING", None, None],
         ["Item", "Choice / Value", "Method / Detail"],
         ["Filtering", ">= (n_rep - 1) values", "synthetic"],
         ["Imputation", "QRILC (MNAR)", "synthetic"],
         ["3. DIFFERENTIAL ANALYSIS (Volcanos)", None, None],
         ["Item", "Choice / Value", "Method / Detail"],
         ["Primary statistic", "limma eBayes", "synthetic"],
         ["Design matrix", "~0 + condition", "synthetic"]])
    keep = set(temps or TEMPS)
    for t, df in _planted(n_proteins, seed).items():
        if t not in keep:
            continue
        if layout == "flat":
            root.mkdir(parents=True, exist_ok=True)
            path = root / f"{t}_ProteomicAnalysis_Results.xlsx"
        else:
            (root / t / "proteomics_output").mkdir(parents=True, exist_ok=True)
            path = root / t / "proteomics_output" / "ProteomicAnalysis_Results.xlsx"
        with pd.ExcelWriter(path, engine="openpyxl") as xw:
            methods.to_excel(xw, sheet_name="Methods", header=False, index=False)
            df.to_excel(xw, sheet_name="Differential_Expression", index=False)
    return root


# ============================================================================= intensities

SIM_TEMPS_C = np.linspace(37, 67, 12)
SIM_CONDS = ["C2", "C1", "Ctl"]


def simulate_intensity_series(n_prot: int = 2000, n_lysates: int = 3, seed: int = 0,
                              frac_shift: float = 0.04, frac_offset: float = 0.02,
                              noise: float = 0.26, lysate_sd: float = 0.25, realistic: bool = True,
                              shift_sizes=(2.0, 4.0, 6.0)):
    """Lysate-treated TPP intensities, one matrix per temperature.

    soluble fraction f(T) = 1 / (1 + exp((T − Tm) / s)), Tm ~ N(52, 6) °C, s ~ U(1.5, 3);
    log2 intensity = base + log2(f + 0.05) + lysate effect + noise. R lysates are split
    into the 3 conditions (Ctl, C1, C2) and the 12 temperatures (paired design).
    C1 stabilises ``frac_shift`` of the proteins (ΔTm in ``shift_sizes``); C2 carries a
    constant +0.8 offset on ``frac_offset`` of the others (handling artefact).

    realistic=True: protein-specific variances (scaled inverse χ², d0 = 2, as limma
    estimates on real TPP runs), t5 noise, ~3 % left-censored values per temperature
    imputed from a narrow normal under the detection limit (QRILC-like).

    Returns (matrices: dict temperature → TemperatureMatrix, truth: DataFrame).
    """
    from .matrices import TemperatureMatrix

    temps = [f"T{i}" for i in range(1, 13)]
    rng = np.random.default_rng(seed)
    names = [f"P{i:05d}" for i in range(n_prot)]
    base = rng.normal(20, 2, n_prot)
    tm = rng.normal(52, 6, n_prot)
    slope = rng.uniform(1.5, 3.0, n_prot)
    if realistic:
        sd_p = np.clip(noise * np.sqrt(2.0 / rng.chisquare(2.0, n_prot)) * 0.6, 0.05, 2.0)
    else:
        sd_p = np.full(n_prot, noise)
    dtm = np.zeros(n_prot)
    shifted = rng.choice(n_prot, int(frac_shift * n_prot), replace=False)
    dtm[shifted] = rng.choice(list(shift_sizes), len(shifted))
    offset = np.zeros(n_prot)
    rest = np.setdiff1d(np.arange(n_prot), shifted)
    offset[rng.choice(rest, int(frac_offset * n_prot), replace=False)] = 0.8
    lys0 = rng.normal(0, lysate_sd, (n_prot, n_lysates))
    mats = {}
    for k, (t, tc) in enumerate(zip(temps, SIM_TEMPS_C)):
        cols, rows = {}, []
        for L in range(n_lysates):
            lys = lys0[:, L] * (1 - 0.04 * k) + rng.normal(0, 0.05, n_prot)
            for cond in SIM_CONDS:
                f = 1.0 / (1.0 + np.exp((tc - (tm + (dtm if cond == "C1" else 0.0))) / slope))
                eps = rng.standard_t(5, n_prot) / np.sqrt(5 / 3) if realistic else rng.normal(0, 1, n_prot)
                y = base + np.log2(f + 0.05) + lys + sd_p * eps
                if cond == "C2":
                    y = y + offset
                lab = f"Sple-{k + 1}{'BCDEFGHIJKLMNOPQ'[L * 3 + SIM_CONDS.index(cond)]}"
                cols[lab] = y
                rows.append({"label": lab, "condition": cond, "replicate": str(L + 1)})
        imp = pd.DataFrame(cols, index=pd.Index(names, name="protein"))
        obs = imp.copy()
        if realistic:
            lod = np.quantile(imp.to_numpy(), 0.03)
            low = imp < lod
            obs = imp.mask(low)
            imp = imp.mask(low, rng.normal(lod - 0.8, 0.3, imp.shape))
        mats[t] = TemperatureMatrix(t, imp, obs, pd.DataFrame(rows), Path("."))
    truth = pd.DataFrame({"name": names, "dTm": dtm, "offset": offset})
    return mats, truth


def write_intensity_runs(root: str | Path, matrices: dict, limma=None) -> Path:
    """Write a Deimos-like layout for the simulated series:
    ``<root>/<T>/proteomics_output/{ProteomicAnalysis_Results.xlsx, matrices/}``.
    The Differential_Expression sheet is computed with Deimos' limma (unpaired,
    as Deimos does by default) when ``limma`` is given, else with Welch t-tests."""
    from scipy import stats as sstats

    root = Path(root)
    contrasts = ["C2_vs_C1", "C2_vs_Ctl", "C1_vs_Ctl"]
    for t, tm in matrices.items():
        out = root / t / "proteomics_output"
        (out / "matrices").mkdir(parents=True, exist_ok=True)
        imp, obs, smp = tm.imputed, tm.observed, tm.samples
        imp.round(6).reset_index().to_csv(out / "matrices" / "log2_imputed.tsv", sep="\t", index=False)
        (obs if obs is not None else imp).round(6).reset_index().to_csv(
            out / "matrices" / "log2_filtered_NA.tsv", sep="\t", index=False)
        smp.to_csv(out / "matrices" / "samples.tsv", sep="\t", index=False)
        de = pd.DataFrame({"name": imp.index, "Protein.Group": imp.index,
                           "Genes": [f"g{n[1:]}" for n in imp.index],
                           "First.Protein.Description": "simulated",
                           "imputed": (obs.isna().any(axis=1) if obs is not None else False),
                           "num_NAs": (obs.isna().sum(axis=1) if obs is not None else 0)})
        if limma is not None:
            from .engine import fit_temperature
            f = fit_temperature(limma, t, imp, smp["condition"].tolist(), smp["replicate"].tolist(),
                                contrasts, paired=False)
            for j, c in enumerate(contrasts):
                de[f"{c}_diff"] = f.diff[:, j]
                de[f"{c}_p.val"] = f.p[:, j]
                de[f"{c}_p.adj"] = np.minimum(1, f.p[:, j] * len(imp))
        else:
            for c in contrasts:
                a, b = c.split("_vs_")
                xa = imp[smp.loc[smp.condition == a, "label"]].to_numpy()
                xb = imp[smp.loc[smp.condition == b, "label"]].to_numpy()
                tt = sstats.ttest_ind(xa, xb, axis=1, equal_var=False)
                de[f"{c}_diff"] = xa.mean(1) - xb.mean(1)
                de[f"{c}_p.val"] = tt.pvalue
                de[f"{c}_p.adj"] = np.minimum(1, tt.pvalue * len(imp))
        de["imputed"] = np.where(de["imputed"], "VRAI", "FAUX")
        with pd.ExcelWriter(out / "ProteomicAnalysis_Results.xlsx", engine="openpyxl") as xw:
            pd.DataFrame([["2. DATA PREPROCESSING", None, None],
                          ["Item", "Choice / Value", "Method / Detail"],
                          ["Imputation", "QRILC (MNAR)", "simulated"],
                          ["Design matrix", "~0 + condition", "simulated"]]).to_excel(
                xw, sheet_name="Methods", header=False, index=False)
            de.to_excel(xw, sheet_name="Differential_Expression", index=False)
    return root
