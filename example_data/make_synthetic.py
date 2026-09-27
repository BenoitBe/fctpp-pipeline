"""Build a synthetic example next to this script.

    python example_data/make_synthetic.py
    fctpp analyze --deimos-root example_data/deimos_runs \
                  --celsius example_data/temperatures.csv -o example_output

Writes 12 Deimos-like workbooks (Differential_Expression + Methods sheets) with
planted thermal shifts (see fctpp/synthetic.py for the list), a temperature
mapping, and raw DIA-NN-like matrices + a template design to try `fctpp prepare`.
"""

from pathlib import Path

import pandas as pd

from fctpp.synthetic import CELSIUS, make_deimos_like_runs, make_raw_pg_matrices

here = Path(__file__).resolve().parent
make_deimos_like_runs(here / "deimos_runs")
pd.DataFrame({"temperature": list(CELSIUS), "celsius": list(CELSIUS.values())}).to_csv(
    here / "temperatures.csv", index=False)
make_raw_pg_matrices(here / "raw", temps=list(CELSIUS))
print(f"Synthetic data written in {here}")
