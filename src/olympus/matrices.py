"""Per-temperature intensity matrices exported by Deimos.

Each Deimos run writes ``<out_dir>/matrices/`` next to its workbook:

* ``log2_imputed.tsv``      proteins × samples, log2, imputed (what Deimos tested)
* ``log2_filtered_NA.tsv``  same proteins, NaN where the value was missing
* ``samples.tsv``           ``label``, ``condition``, ``replicate`` (+ ``subject`` if any)

These are what the shift test refits. The ``protein`` key is Deimos' ``name``
(first accession of Protein.Group), the same key as the Differential_Expression
sheet.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)


@dataclass
class TemperatureMatrix:
    temperature: str
    imputed: pd.DataFrame            # proteins × samples (log2, imputed)
    observed: pd.DataFrame | None    # same, NaN where missing
    samples: pd.DataFrame            # label, condition, replicate[, subject]
    source: Path

    @property
    def units(self) -> list[str]:
        """'condition|block' per sample column, block = subject or replicate."""
        blk = self.block_labels()
        return [f"{c}|{b}" for c, b in zip(self.samples["condition"], blk)]

    def block_labels(self) -> list[str]:
        s = self.samples
        if "subject" in s.columns and s["subject"].notna().all() and (s["subject"].astype(str) != "").all():
            return s["subject"].astype(str).tolist()
        return s["replicate"].astype(str).tolist()


def matrices_dir(results_path: str | Path) -> Path:
    return Path(results_path).parent / "matrices"


def has_matrices(results_path: str | Path) -> bool:
    d = matrices_dir(results_path)
    return (d / "log2_imputed.tsv").exists() and (d / "samples.tsv").exists()


def _read_matrix(path: Path, labels: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    key = "protein" if "protein" in df.columns else ("name" if "name" in df.columns else df.columns[0])
    df[key] = df[key].astype(str).str.strip()
    missing = [lab for lab in labels if lab not in df.columns]
    if missing:
        raise ValueError(f"{path}: sample column(s) {missing} listed in samples.tsv are absent")
    m = df.set_index(key)[labels].apply(pd.to_numeric, errors="coerce")
    if m.index.duplicated().any():
        m = m[~m.index.duplicated(keep="first")]
    return m


def load_temperature_matrix(temperature: str, results_path: str | Path) -> TemperatureMatrix:
    d = matrices_dir(results_path)
    if not has_matrices(results_path):
        raise FileNotFoundError(
            f"[{temperature}] no Deimos matrices in {d} (log2_imputed.tsv + samples.tsv). "
            f"Rerun Deimos with export_matrices: true (default since Deimos e160463).")
    samples = pd.read_csv(d / "samples.tsv", sep="\t", dtype=str)
    for col in ("label", "condition", "replicate"):
        if col not in samples.columns:
            raise ValueError(f"{d / 'samples.tsv'}: column '{col}' missing")
    samples = samples.astype(object).where(samples.notna(), None)
    labels = samples["label"].astype(str).tolist()
    imputed = _read_matrix(d / "log2_imputed.tsv", labels)
    observed = None
    if (d / "log2_filtered_NA.tsv").exists():
        observed = _read_matrix(d / "log2_filtered_NA.tsv", labels).reindex(imputed.index)
    if imputed.isna().any().any():
        n = int(imputed.isna().any(axis=1).sum())
        log.warning("[%s] %d protein(s) with missing values in log2_imputed.tsv — dropped from the "
                    "shift test", temperature, n)
        imputed = imputed.dropna()
    return TemperatureMatrix(temperature, imputed, observed, samples, d)


def load_series_matrices(runs: dict[str, str | Path]) -> dict[str, TemperatureMatrix]:
    return {t: load_temperature_matrix(t, p) for t, p in runs.items()}
