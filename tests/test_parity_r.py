"""Parity with the R package on its own outputs (needs the FC.TPP R project).

    FCTPP_R_PROJECT=/path/to/FC.TPP pytest -m realdata

1. Call logic: the R ``TPP_data`` sheet re-analysed with ``rule=r_legacy`` gives
   identical ``_logpval`` and ``_threshold`` columns.
2. Presence statistics: recomputed from the raw ``TPP_data/T*.tsv`` (cRAP
   removal + DEP ``filter_missval(thr = 1)`` + DEP-style label matching), the
   ``nb_temp`` / ``max_consecutive`` values equal the R ``Stats_proteins`` sheet.
"""

import numpy as np
import pandas as pd
import pytest

from fctpp.design import match_samples, read_design
from fctpp.inputs import load_long_table
from fctpp.prepare import discover_raw_files, split_pg_matrix
from fctpp.temperatures import TemperatureSeries
from fctpp.tpp import CallParams, presence_stats, run_tpp

pytestmark = pytest.mark.realdata


def test_call_logic_matches_r(r_project):
    li = load_long_table(r_project / "TPP_results.xlsx", "TPP_data")
    res = run_tpp(li.results, call_params=CallParams(rule="r_legacy"), contrasts=li.contrasts)
    r = pd.read_excel(r_project / "TPP_results.xlsx", sheet_name="TPP_data")
    m = r.merge(res.data, on=["Temperature", "name"], suffixes=("_R", "_py"))
    assert len(m) == len(r) == len(res.data)
    for c in li.contrasts:
        assert (m[f"{c}_threshold_R"] == m[f"{c}_threshold_py"]).all()
        np.testing.assert_allclose(m[f"{c}_logpval_R"], m[f"{c}_logpval_py"], rtol=0, atol=1e-12)


def dep_filter_missval(x: pd.DataFrame, conditions: np.ndarray, thr: int = 1) -> pd.Series:
    """DEP::filter_missval: <= thr missing values in at least one condition."""
    miss = pd.DataFrame({c: x.loc[:, conditions == c].isna().sum(axis=1) for c in np.unique(conditions)})
    return (miss <= thr).any(axis=1)


def test_presence_statistics_match_r(r_project):
    raw = r_project / "TPP_data"
    design = read_design(raw / "experimental_design.xlsx")
    parts = []
    for t, p in discover_raw_files(raw).items():
        df = pd.read_csv(p, sep="\t")
        crap = (df["Protein.Group"].astype(str).str.startswith("cRAP")
                | df["Protein.Names"].astype(str).str.startswith("cRAP"))
        df = df[~crap]
        _, samples = split_pg_matrix(df)
        m = match_samples(samples, design)
        x = df[m["column"]].replace(0, np.nan)
        keep = dep_filter_missval(x, m["condition"].to_numpy())
        names = df.loc[keep, "Protein.Group"].astype(str).str.split(";").str[0]
        parts.append(pd.DataFrame({"Temperature": t, "name": names}))
    long = pd.concat(parts)
    st = presence_stats(long, TemperatureSeries.from_ids(long["Temperature"].unique()))
    r = pd.read_excel(r_project / "TPP_results.xlsx", sheet_name="Stats_proteins")
    m = r.merge(st, on="name", how="outer", suffixes=("_R", "_py"), indicator=True)
    assert (m["_merge"] == "both").all()
    assert (m["nb_temp_R"] == m["nb_temp_py"]).all()
    assert (m["max_consecutive_R"] == m["max_consecutive_py"]).all()
