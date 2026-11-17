"""Parity with the R package on its own outputs (needs the R FC.TPP project).

    OLYMPUS_R_PROJECT=/path/to/FC.TPP pytest -m realdata

1. Call logic: the R ``TPP_data`` sheet re-analysed with ``rule=r_legacy`` gives
   identical ``_logpval`` and ``_threshold`` columns.
2. Presence statistics: recomputed from the raw ``TPP_data/T*.tsv`` with the DEP
   rules (cRAP removal, ``filter_missval(thr = 1)``, label matching after removal
   of the longest common prefix as in ``DEP::make_se``), ``nb_temp`` and
   ``max_consecutive`` equal the R ``Stats_proteins`` sheet.
"""

import os
import re

import numpy as np
import pandas as pd
import pytest

from olympus.inputs import load_long_table
from olympus.temperatures import TemperatureSeries, natural_key
from olympus.tpp import CallParams, presence_stats, run_tpp

pytestmark = pytest.mark.realdata

ANNOTATIONS = {"Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
               "First.Protein.Description", "Proteotypic", "N.Sequences", "N.Proteotypic.Sequences"}


def test_call_logic_matches_r(r_project):
    li = load_long_table(r_project / "TPP_results.xlsx", "TPP_data")
    res = run_tpp(li.results, call_params=CallParams(rule="r_legacy"), contrasts=li.contrasts)
    r = pd.read_excel(r_project / "TPP_results.xlsx", sheet_name="TPP_data")
    m = r.merge(res.data, on=["Temperature", "name"], suffixes=("_R", "_py"))
    assert len(m) == len(r) == len(res.data)
    for c in li.contrasts:
        assert (m[f"{c}_threshold_R"] == m[f"{c}_threshold_py"]).all()
        np.testing.assert_allclose(m[f"{c}_logpval_R"], m[f"{c}_logpval_py"], rtol=0, atol=1e-12)


def _suffixes(labels):
    """DEP delete_prefix + make.names, as in make_se()."""
    canon = [re.sub(r"[^A-Za-z0-9._]", ".", str(x)) for x in labels]
    prefix = os.path.commonprefix(canon)
    return [x[len(prefix):] for x in canon]


def test_presence_statistics_match_r(r_project):
    raw = r_project / "TPP_data"
    design = pd.read_excel(raw / "experimental_design.xlsx", dtype=str)
    d_suffix = dict(zip(_suffixes(design["label"]), design["condition"]))
    parts = []
    for p in sorted(raw.glob("T*.tsv"), key=lambda f: natural_key(f.stem)):
        df = pd.read_csv(p, sep="\t")
        crap = (df["Protein.Group"].astype(str).str.startswith("cRAP")
                | df["Protein.Names"].astype(str).str.startswith("cRAP"))
        df = df[~crap]
        samples = [c for c in df.columns if c not in ANNOTATIONS]
        clean = [re.search(r"(Sple-[^_]+)_", c).group(1) for c in samples]
        cond = np.array([d_suffix[s] for s in _suffixes(clean)])
        x = df[samples].replace(0, np.nan)
        miss = pd.DataFrame({c: x.loc[:, cond == c].isna().sum(axis=1) for c in np.unique(cond)})
        keep = (miss <= 1).any(axis=1)           # DEP::filter_missval(thr = 1)
        names = df.loc[keep, "Protein.Group"].astype(str).str.split(";").str[0]
        parts.append(pd.DataFrame({"Temperature": p.stem, "name": names}))
    long = pd.concat(parts)
    st = presence_stats(long, TemperatureSeries.from_ids(long["Temperature"].unique()))
    r = pd.read_excel(r_project / "TPP_results.xlsx", sheet_name="Stats_proteins")
    m = r.merge(st, on="name", how="outer", suffixes=("_R", "_py"), indicator=True)
    assert (m["_merge"] == "both").all()
    assert (m["nb_temp_R"] == m["nb_temp_py"]).all()
    assert (m["max_consecutive_R"] == m["max_consecutive_py"]).all()
