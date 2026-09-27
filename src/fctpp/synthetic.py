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


def make_deimos_like_runs(root: str | Path, n_proteins: int = 300, seed: int = 7) -> Path:
    """``<root>/<T>/deimos_output/ProteomicAnalysis_Results.xlsx`` for T1..T12."""
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
    for t, df in _planted(n_proteins, seed).items():
        out = root / t / "deimos_output"
        out.mkdir(parents=True, exist_ok=True)
        with pd.ExcelWriter(out / "ProteomicAnalysis_Results.xlsx", engine="openpyxl") as xw:
            methods.to_excel(xw, sheet_name="Methods", header=False, index=False)
            df.to_excel(xw, sheet_name="Differential_Expression", index=False)
    return root


def make_raw_pg_matrices(root: str | Path, n_proteins: int = 60, temps: list[str] | None = None,
                         seed: int = 3, four_annotation_columns: bool = True) -> tuple[Path, Path]:
    """Raw DIA-NN-like pg_matrix per temperature + R-style template design.

    Sample columns follow the TPP run names: ``...\\250429-Sple-<n><L>_Slot1-..d``
    with letters B..J = Ctl/C1/C2 × 3 replicates as in the FC.TPP design.
    """
    rng = np.random.default_rng(seed)
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    temps = temps or TEMPS[:4]
    letters = list("BCDEFGHIJ")
    cond = {"B": "Ctl", "C": "C1", "D": "C2", "E": "Ctl", "F": "C1", "G": "C2",
            "H": "Ctl", "I": "C1", "J": "C2"}
    rep = {L: str(i // 3 + 1) for i, L in enumerate(letters)}
    base = rng.normal(22, 2, n_proteins)
    for t in temps:
        n = int(t[1:])
        cols = {}
        for k, L in enumerate(letters):
            x = 2 ** (base + rng.normal(0, 0.2, n_proteins) - 0.3 * n)
            x[rng.random(n_proteins) < 0.05] = np.nan
            cols[f"E:\\DATA\\TPP\\250429-Sple-{n}{L}_Slot1-{k + 10}_1_{300 + k}.d"] = x
        ann = {"Protein.Group": [f"Q{i:05d}" for i in range(n_proteins)],
               "Protein.Names": [f"PRT{i}_STAA8" for i in range(n_proteins)],
               "Genes": [f"g{i}" for i in range(n_proteins)],
               "First.Protein.Description": ["synthetic"] * n_proteins}
        ann["Protein.Group"][0] = "cRAP-KERATIN"
        if not four_annotation_columns:
            ann["N.Sequences"] = rng.integers(1, 20, n_proteins)
            ann["N.Proteotypic.Sequences"] = rng.integers(1, 20, n_proteins)
        pd.DataFrame({**ann, **cols}).to_csv(root / f"{t}.tsv", sep="\t", index=False)
    design = pd.DataFrame({"label": [f"Sple.1{L}" for L in "DCBGFEJIH"],
                           "condition": [cond[L] for L in "DCBGFEJIH"],
                           "replicate": [rep[L] for L in "DCBGFEJIH"]})
    design_path = root / "experimental_design.csv"
    design.to_csv(design_path, sep=";", index=False)
    return root, design_path
