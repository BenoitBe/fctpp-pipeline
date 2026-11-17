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


@pytest.fixture(scope="session")
def limma_engine():
    """Deimos' limma_ebayes (DEIMOS_DIR); the shift-test suites are skipped without it."""
    from olympus.engine import load_limma
    try:
        return load_limma(os.environ.get("DEIMOS_DIR"))
    except ImportError:
        pytest.skip("set DEIMOS_DIR to a Deimos clone (limma_ebayes.py) to run the shift-test suites")


@pytest.fixture(scope="session")
def sim_series():
    from olympus.synthetic import simulate_intensity_series
    return simulate_intensity_series(n_prot=600, n_lysates=3, seed=1)
