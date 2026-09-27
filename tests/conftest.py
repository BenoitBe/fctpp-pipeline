import os
from pathlib import Path

import pytest

from fctpp.synthetic import make_deimos_like_runs, make_raw_pg_matrices


@pytest.fixture(scope="session")
def deimos_root(tmp_path_factory) -> Path:
    return make_deimos_like_runs(tmp_path_factory.mktemp("deimos_runs"))


@pytest.fixture()
def raw_dir(tmp_path) -> tuple[Path, Path]:
    return make_raw_pg_matrices(tmp_path / "raw")


@pytest.fixture(scope="session")
def r_project() -> Path:
    p = os.environ.get("FCTPP_R_PROJECT")
    if not p or not Path(p, "TPP_results.xlsx").exists():
        pytest.skip("set FCTPP_R_PROJECT to the FC.TPP R project folder (TPP_results.xlsx, TPP_data/)")
    return Path(p)


@pytest.fixture(scope="session")
def deimos_dir() -> Path:
    p = os.environ.get("FCTPP_DEIMOS_DIR")
    if not p or not Path(p, "deimos.py").exists():
        pytest.skip("set FCTPP_DEIMOS_DIR to a clone of deimos-proteomics-pipeline")
    return Path(p)
