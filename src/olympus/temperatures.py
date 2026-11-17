"""Temperature series: ordering, adjacency and axis labels.

In the R package the series was hard-coded as ``T1..T12`` (factor levels in
``tidy_tpp`` and the 3D threshold planes in ``plot_and_3D``). Here the series is
whatever the inputs contain, ordered naturally (T1 < T2 < ... < T10) or by the
temperature in °C when a mapping is supplied. "Consecutive" always means
adjacent in this ordered series.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping

import pandas as pd

_SPLIT = re.compile(r"(\d+(?:\.\d+)?)")


def natural_key(value: object) -> list:
    """Sort key so that 'T2' < 'T10' (equivalent to stringr::str_sort(numeric=TRUE))."""
    parts = _SPLIT.split(str(value))
    return [float(p) if i % 2 else p.lower() for i, p in enumerate(parts)]


def natural_sorted(values: Iterable[object]) -> list[str]:
    return sorted((str(v) for v in values), key=natural_key)


@dataclass(frozen=True)
class TemperatureSeries:
    """Ordered temperature identifiers with optional values in °C."""

    ids: tuple[str, ...]
    celsius: Mapping[str, float] | None = None
    _pos: dict = field(default_factory=dict, repr=False, compare=False)

    def __post_init__(self) -> None:
        if len(set(self.ids)) != len(self.ids):
            raise ValueError(f"Duplicated temperature identifiers: {self.ids}")
        object.__setattr__(self, "_pos", {t: i for i, t in enumerate(self.ids)})

    # ------------------------------------------------------------------ build
    @classmethod
    def from_ids(cls, ids: Iterable[object],
                 celsius: Mapping[str, float] | None = None) -> "TemperatureSeries":
        uniq = list(dict.fromkeys(str(i) for i in ids))
        if celsius:
            cel = {str(k): float(v) for k, v in celsius.items()}
            missing = [t for t in uniq if t not in cel]
            if missing:
                raise ValueError(
                    f"No temperature in °C for {missing}. Either give a value for "
                    f"every temperature or none at all.")
            ordered = sorted(uniq, key=lambda t: (cel[t], natural_key(t)))
            values = [cel[t] for t in ordered]
            if len(set(values)) != len(values):
                raise ValueError(f"Two temperatures share the same °C value: {cel}")
            return cls(tuple(ordered), {t: cel[t] for t in ordered})
        return cls(tuple(natural_sorted(uniq)), None)

    # ---------------------------------------------------------------- queries
    def __len__(self) -> int:
        return len(self.ids)

    def __iter__(self):
        return iter(self.ids)

    def position(self, t: object) -> int:
        return self._pos[str(t)]

    def positions(self, values: Iterable[object]) -> list[int]:
        return [self._pos[str(v)] for v in values]

    @property
    def reference(self) -> str:
        """Lowest temperature of the series (abundance reference)."""
        return self.ids[0]

    def x(self, t: object) -> float:
        """Numeric x coordinate: °C if known, else 1-based rank."""
        t = str(t)
        if self.celsius:
            return float(self.celsius[t])
        return float(self._pos[t] + 1)

    def label(self, t: object) -> str:
        t = str(t)
        if self.celsius:
            return f"{t} ({self.celsius[t]:g} °C)"
        return t

    @property
    def axis_title(self) -> str:
        return "Temperature (°C)" if self.celsius else "Temperature"

    def tick_labels(self) -> list[str]:
        if self.celsius:
            return [f"{self.celsius[t]:g}" for t in self.ids]
        return list(self.ids)

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame({
            "Temperature": list(self.ids),
            "order": range(1, len(self.ids) + 1),
            "Temperature_C": [self.celsius[t] if self.celsius else None for t in self.ids],
        })


def read_celsius_mapping(source: str | Path | Mapping | None) -> dict[str, float] | None:
    """Temperature → °C mapping from a dict (YAML) or a CSV/TSV/XLSX file.

    File columns: ``temperature`` (identifier, e.g. T1) and ``celsius``.
    """
    if source is None or source == "":
        return None
    if isinstance(source, Mapping):
        return {str(k): float(v) for k, v in source.items()} or None
    path = Path(source)
    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path, sep=None, engine="python")
    cols = {c.lower().strip(): c for c in df.columns}
    t_col = cols.get("temperature") or cols.get("temperature_id") or cols.get("id")
    c_col = cols.get("celsius") or cols.get("temperature_c") or cols.get("degc")
    if t_col is None or c_col is None:
        raise ValueError(f"{path}: expected columns 'temperature' and 'celsius', got {list(df.columns)}")
    return {str(t).strip(): float(c) for t, c in zip(df[t_col], df[c_col])}
