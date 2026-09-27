"""End-to-end with the real Deimos (slow, ~30 s per temperature).

    FCTPP_DEIMOS_DIR=../deimos-proteomics-pipeline FCTPP_R_PROJECT=../FC.TPP pytest -m deimos

Runs prepare → Deimos → analyze on three temperatures of the real TPP data and
checks that, for proteins without imputed values, Deimos' fold changes equal
those of the R package (identical log2 intensities, identical ~0 + condition
limma model).
"""

import shutil

import numpy as np
import pandas as pd
import pytest

from fctpp.cli import main

pytestmark = pytest.mark.deimos


def test_prepare_run_analyze_on_real_data(tmp_path, deimos_dir, r_project):
    raw = tmp_path / "raw"
    raw.mkdir()
    for t in ("T1", "T2", "T3"):
        shutil.copy(r_project / "TPP_data" / f"{t}.tsv", raw)
    design = r_project / "TPP_data" / "experimental_design.xlsx"
    runs, out = tmp_path / "runs", tmp_path / "out"
    code = main(["all", "--raw-dir", str(raw), "--design", str(design), "--runs", str(runs),
                 "--deimos-dir", str(deimos_dir), "-o", str(out), "--no-plots",
                 "--min-temperatures", "3", "--min-streak", "3",
                 "--set", "prepare.deimos.control_condition=Ctl"])
    assert code == 0
    data = pd.read_csv(out / "TPP_data.csv")
    assert sorted(data["Temperature"].unique()) == ["T1", "T2", "T3"]
    assert set(pd.read_csv(runs / "runs_status.csv")["status"]) == {"ok"}

    r = pd.read_excel(r_project / "TPP_results.xlsx", sheet_name="TPP_data")
    m = r.merge(data, on=["Temperature", "name"], suffixes=("_R", "_py"))
    complete = m[m["num_NAs"] == 0]
    assert len(complete) > 5000
    for c in ("C1_vs_Ctl", "C2_vs_C1", "C2_vs_Ctl"):
        np.testing.assert_allclose(complete[f"{c}_diff_R"], complete[f"{c}_diff_py"], atol=1e-10)
