"""FC.TPP — fold-change based Thermal Proteome Profiling (Python port).

Per-temperature differential statistics are produced upstream by Deimos
(one Deimos run per temperature). FC.TPP assembles those runs, filters proteins
on their temperature coverage and calls thermal-shift candidates from the
fold-change / p-value profiles across the temperature series.
"""

__version__ = "0.1.0"

from .temperatures import TemperatureSeries  # noqa: E402,F401
from .tpp import CallParams, FilterParams, TPPResult, run_tpp  # noqa: E402,F401
