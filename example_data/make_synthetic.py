"""Build a synthetic example next to this script.

    python example_data/make_synthetic.py
    olympus analyze --deimos-root example_data/deimos_results \
                  --celsius example_data/temperatures.csv -o example_output

Writes 12 Deimos-like result workbooks (Differential_Expression + Methods
sheets), one per temperature, with planted thermal shifts (see
olympus/synthetic.py for the list), and a temperature → °C mapping.
"""

from pathlib import Path

import pandas as pd

from olympus.synthetic import CELSIUS, make_deimos_like_runs

here = Path(__file__).resolve().parent
make_deimos_like_runs(here / "deimos_results")
pd.DataFrame({"temperature": list(CELSIUS), "celsius": list(CELSIUS.values())}).to_csv(
    here / "temperatures.csv", index=False)
print(f"Synthetic data written in {here}")
