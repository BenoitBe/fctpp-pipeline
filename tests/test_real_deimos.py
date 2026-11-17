"""Real Deimos result workbooks vs the R package (slow; needs both folders).

    OLYMPUS_DEIMOS_RESULTS=/path/deimos_results OLYMPUS_R_PROJECT=/path/FC.TPP pytest -m deimos

For proteins without imputed values, the Deimos log2FC must equal those of the
R/DEP pipeline (same log2 intensities, same ~0 + condition limma model), and the
two R hits must be found again.
"""

import numpy as np
import pandas as pd
import pytest

from olympus.cli import main

pytestmark = pytest.mark.deimos


def test_real_deimos_results_match_r(tmp_path, real_deimos_results, r_project):
    out = tmp_path / "out"
    assert main(["analyze", "--deimos-root", str(real_deimos_results), "-o", str(out), "--no-plots"]) == 0
    data = pd.read_csv(out / "TPP_data.csv")
    r = pd.read_excel(r_project / "TPP_results.xlsx", sheet_name="TPP_data")
    m = r.merge(data, on=["Temperature", "name"], suffixes=("_R", "_py"))
    complete = m[m["num_NAs"] == 0]
    assert len(complete) > 0.9 * len(m)
    for c in ("C1_vs_Ctl", "C2_vs_C1", "C2_vs_Ctl"):
        np.testing.assert_allclose(complete[f"{c}_diff_R"], complete[f"{c}_diff_py"], atol=1e-10)
    hits = pd.read_csv(out / "TPP_hits.csv")
    r_hits = set(r.loc[r["C1_vs_Ctl_threshold"] == "Over", "name"])
    assert r_hits <= set(hits.loc[hits["contrast"] == "C1_vs_Ctl", "name"])
