import os
from pathlib import Path

import pytest

from olympus.synthetic import make_deimos_like_runs


@pytest.fixture(scope="session")
def deimos_root(tmp_path_factory) -> Path:
    """12 synthetic Deimos result workbooks (<root>/<T>/proteomics_output/...)."""
    return make_deimos_like_runs(tmp_path_factory.mktemp("deimos_results"))


@pytest.fixture(scope="session")
def r_project() -> Path:
    p = os.environ.get("OLYMPUS_R_PROJECT")
    if not p or not Path(p, "TPP_results.xlsx").exists():
        pytest.skip("set OLYMPUS_R_PROJECT to the R FC.TPP project folder (TPP_results.xlsx, TPP_data/)")
    return Path(p)


@pytest.fixture(scope="session")
def real_deimos_results() -> Path:
    p = os.environ.get("OLYMPUS_DEIMOS_RESULTS")
    if not p or not Path(p).is_dir():
        pytest.skip("set OLYMPUS_DEIMOS_RESULTS to a folder holding the real Deimos result workbooks")
    return Path(p)
