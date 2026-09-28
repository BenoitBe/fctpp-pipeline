"""Readers: Deimos per-temperature workbooks, or one long table.

Deimos writes ``ProteomicAnalysis_Results.xlsx`` with a ``Differential_Expression``
sheet: ``name``, ``Protein.Group``, ``Genes``, ``First.Protein.Description``,
``imputed``, ``num_NAs`` and, per contrast, ``<c>_diff``, ``<c>_p.val``,
``<c>_p.adj`` (+ ``Pi_Score_<c>``, ``Robustness_Score_<c>``). ``name`` is the
first accession of Protein.Group — the same key as ``name`` in the R package.
"""

from __future__ import annotations

import hashlib
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .temperatures import natural_key
from .tpp import detect_contrasts

log = logging.getLogger(__name__)

RESULTS_WORKBOOK = "ProteomicAnalysis_Results.xlsx"
DE_SHEET = "Differential_Expression"
ANNOTATIONS = ("Protein.Group", "Genes", "First.Protein.Description", "N.Sequences",
               "imputed", "num_NAs")
# Deimos 'Methods' items compared across temperatures (processing must be identical)
METHOD_KEYS = ("Filtering", "Imputation", "Primary statistic", "Design matrix",
               "Contrasts", "Multiple correction", "Proteins retained")


@dataclass
class LoadedInputs:
    results: pd.DataFrame                   # long table, one row per (Temperature, name)
    contrasts: list[str]
    provenance: pd.DataFrame                # one row per temperature
    methods: pd.DataFrame = field(default_factory=pd.DataFrame)   # Deimos Methods, per temperature
    warnings: list[str] = field(default_factory=list)


def sha256(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _as_bool(v) -> object:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    s = str(v).strip().lower()
    if s in ("true", "vrai", "1", "yes", "oui", "t"):
        return True
    if s in ("false", "faux", "0", "no", "non", "f"):
        return False
    return None


# ----------------------------------------------------------------------------
# Discovery
# ----------------------------------------------------------------------------

DEFAULT_WORKBOOK_GLOB = "**/*ProteomicAnalysis_Results*.xlsx"


def temperature_id(rel_path: Path, temperature_regex: str | None = None) -> str | None:
    """Temperature identifier of a workbook, from its path relative to the root.

    * with ``temperature_regex``: first capture group (or whole match) searched in
      the relative path, e.g. ``(T\\d+)`` on ``TPP_T3_output/ProteomicAnalysis_Results.xlsx``
      gives ``T3``; paths that do not match are skipped (returns None);
    * otherwise: the first folder under the root (``<root>/T3/.../file.xlsx`` → ``T3``),
      or the file name without the Deimos suffix for a flat folder
      (``T3_ProteomicAnalysis_Results.xlsx`` → ``T3``).
    """
    rel = rel_path.as_posix()
    if temperature_regex:
        m = re.search(temperature_regex, rel)
        if not m:
            return None
        return m.group(1) if m.groups() else m.group(0)
    if len(rel_path.parts) > 1:
        return rel_path.parts[0]
    stem = re.sub(r"[_\-. ]*ProteomicAnalysis_Results[_\-. ]*", "", rel_path.stem)
    return stem or None


def discover_deimos_runs(root: str | Path, temperature_regex: str | None = None,
                         pattern: str = DEFAULT_WORKBOOK_GLOB) -> dict[str, Path]:
    """Find the X Deimos result workbooks (one per temperature) under ``root``.

    Accepted layouts (any number of temperatures)::

        root/T1/ProteomicAnalysis_Results.xlsx            (Deimos out_dir per temperature)
        root/T1/proteomics_output/ProteomicAnalysis_Results.xlsx
        root/T1_ProteomicAnalysis_Results.xlsx            (renamed copies in one folder)

    Excel lock files (``~$...``) are ignored. Two workbooks resolving to the same
    temperature (e.g. a Deimos per-organism split) stop the run: use
    ``temperature_regex`` or a manifest to disambiguate.
    """
    root = Path(root)
    if not root.is_dir():
        raise FileNotFoundError(f"Deimos results folder not found: {root}")
    found: dict[str, Path] = {}
    for wb in sorted(root.glob(pattern)):
        if wb.name.startswith("~$") or not wb.is_file():
            continue
        temp = temperature_id(wb.relative_to(root), temperature_regex)
        if temp is None:
            continue
        if temp in found:
            raise ValueError(f"Two Deimos workbooks for temperature '{temp}': {found[temp]} and {wb} "
                             f"(set input.temperature_regex or use a manifest)")
        found[temp] = wb
    if not found:
        raise FileNotFoundError(f"No Deimos workbook matching '{pattern}' under {root}")
    return dict(sorted(found.items(), key=lambda kv: natural_key(kv[0])))


def read_manifest(path: str | Path) -> tuple[dict[str, Path], dict[str, float] | None]:
    """CSV listing the X Deimos workbooks: columns ``temperature``,
    ``results_path`` and optionally ``temperature_c`` (relative paths are
    resolved from the manifest folder)."""
    path = Path(path)
    df = pd.read_csv(path, sep=None, engine="python")
    cols = {c.lower(): c for c in df.columns}
    if "temperature" not in cols or "results_path" not in cols:
        raise ValueError(f"{path}: expected columns 'temperature' and 'results_path'")
    base = path.parent
    runs = {}
    for t, p in zip(df[cols["temperature"]], df[cols["results_path"]]):
        p = Path(str(p))
        runs[str(t)] = p if p.is_absolute() else (base / p)
    celsius = None
    if "temperature_c" in cols and df[cols["temperature_c"]].notna().all():
        celsius = {str(t): float(c) for t, c in zip(df[cols["temperature"]], df[cols["temperature_c"]])}
    return runs, celsius


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------

def read_deimos_methods(path: str | Path) -> dict[str, str]:
    """``Item → Choice / Value`` from the Deimos 'Methods' sheet (section-prefixed
    when an item appears in several sections)."""
    try:
        m = pd.read_excel(path, sheet_name="Methods", header=None, dtype=str)
    except (ValueError, KeyError):
        return {}
    out: dict[str, str] = {}
    section = ""
    for _, row in m.iterrows():
        a = row.iloc[0] if len(row) > 0 else None
        b = row.iloc[1] if len(row) > 1 else None
        if isinstance(a, str) and re.match(r"^\d+\.\s", a) and (b is None or pd.isna(b)):
            section = a.strip()
            continue
        if isinstance(a, str) and isinstance(b, str) and a != "Item":
            key = a.strip()
            if key in out:
                key = f"{section} / {key}"
            out[key] = b.strip()
    return out


def read_deimos_workbook(path: str | Path, sheet: str = DE_SHEET) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Deimos workbook not found: {path}")
    try:
        df = pd.read_excel(path, sheet_name=sheet)
    except ValueError as e:
        raise ValueError(f"{path}: no '{sheet}' sheet — is it a Deimos output?") from e
    if "name" not in df.columns:
        raise ValueError(f"{path}: '{sheet}' has no 'name' column")
    return df


def standardise_table(df: pd.DataFrame, temperature: str | None = None,
                      contrasts: list[str] | None = None) -> pd.DataFrame:
    """Keep key, annotations and ``<c>_{diff,p.val,p.adj}`` (+ robustness)."""
    contrasts = contrasts or detect_contrasts(df.columns)
    keep = ["name"] + [a for a in ANNOTATIONS if a in df.columns]
    for c in contrasts:
        for s in ("diff", "p.val", "p.adj"):
            if f"{c}_{s}" in df.columns:
                keep.append(f"{c}_{s}")
        if f"Robustness_Score_{c}" in df.columns:
            keep.append(f"Robustness_Score_{c}")
    out = df[keep].copy()
    out["name"] = out["name"].astype(str).str.strip()
    if "imputed" in out.columns:
        out["imputed"] = out["imputed"].map(_as_bool)
    if temperature is not None:
        out.insert(0, "Temperature", str(temperature))
    return out


def harmonise_contrasts(tables: dict[str, pd.DataFrame]) -> tuple[dict[str, pd.DataFrame], list[str], list[str]]:
    """Make every temperature expose the same contrast names.

    If the reference name ``A_vs_B`` is missing at one temperature but ``B_vs_A``
    exists (e.g. the Deimos control auto-detection behaved differently), the
    contrast is flipped (diff × -1, p-values unchanged). Contrasts that cannot be
    reconciled are dropped with a warning.
    """
    warnings: list[str] = []
    per_temp = {t: detect_contrasts(df.columns) for t, df in tables.items()}
    first = next(iter(per_temp))
    reference = list(per_temp[first])
    out: dict[str, pd.DataFrame] = {}
    final = list(reference)
    for t, df in tables.items():
        df = df.copy()
        have = set(per_temp[t])
        used = set()
        for c in reference:
            if c in have:
                used.add(c)
                continue
            if "_vs_" in c:
                a, b = c.split("_vs_", 1)
                flipped = f"{b}_vs_{a}"
                if flipped in have:
                    used.add(flipped)
                    df[f"{c}_diff"] = -pd.to_numeric(df[f"{flipped}_diff"], errors="coerce")
                    for s in ("p.val", "p.adj"):
                        if f"{flipped}_{s}" in df.columns:
                            df[f"{c}_{s}"] = df[f"{flipped}_{s}"]
                    if f"Robustness_Score_{flipped}" in df.columns:
                        df[f"Robustness_Score_{c}"] = df[f"Robustness_Score_{flipped}"]
                    warnings.append(f"[{t}] contrast '{flipped}' flipped to '{c}' (diff × -1)")
                    continue
            warnings.append(f"[{t}] contrast '{c}' missing — dropped from the whole analysis")
            if c in final:
                final.remove(c)
        extra = have - used
        if extra:
            warnings.append(f"[{t}] contrast(s) {sorted(extra)} absent at {first} — ignored")
        out[t] = df
    return out, final, warnings


def load_deimos_runs(runs: dict[str, str | Path]) -> LoadedInputs:
    """Read one Deimos workbook per temperature into a long table."""
    tables, prov, methods = {}, [], {}
    for t, p in runs.items():
        p = Path(p)
        df = read_deimos_workbook(p)
        tables[t] = df
        m = read_deimos_methods(p)
        methods[t] = m
        prov.append({"Temperature": t, "path": str(p.resolve()), "sha256": sha256(p),
                     "n_proteins": len(df), "contrasts": ", ".join(detect_contrasts(df.columns)),
                     **{f"deimos:{k}": m.get(k) for k in METHOD_KEYS}})
        log.info("[%s] %d proteins read from %s", t, len(df), p)

    tables, contrasts, warns = harmonise_contrasts(tables)
    if not contrasts:
        raise ValueError("No contrast shared by all temperatures")
    long = pd.concat([standardise_table(df, t, contrasts) for t, df in tables.items()],
                     ignore_index=True)

    methods_df = pd.DataFrame(methods).rename_axis("item").reset_index()
    for key in METHOD_KEYS:
        vals = {t: m.get(key) for t, m in methods.items() if m.get(key) is not None}
        distinct = set(vals.values())
        if key != "Proteins retained" and len(distinct) > 1:
            warns.append(f"Deimos '{key}' differs between temperatures: {vals} — "
                         f"the temperatures were not processed identically")
    for w in warns:
        log.warning(w)
    return LoadedInputs(long, contrasts, pd.DataFrame(prov), methods_df, warns)


def load_long_table(path: str | Path, sheet: str | int | None = None) -> LoadedInputs:
    """One table already in long format (``Temperature``, ``name``, ``<c>_diff``,
    ``<c>_p.val``...), e.g. the ``TPP_data`` sheet of the R package."""
    path = Path(path)
    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path, sheet_name=sheet if sheet is not None else 0)
    else:
        df = pd.read_csv(path, sep=None, engine="python")
    if "Temperature" not in df.columns:
        raise ValueError(f"{path}: long table needs a 'Temperature' column")
    contrasts = detect_contrasts(df.columns)
    parts = [standardise_table(sub, str(t), contrasts) for t, sub in df.groupby("Temperature", sort=False)]
    long = pd.concat(parts, ignore_index=True)
    prov = pd.DataFrame([{"Temperature": "all", "path": str(path.resolve()), "sha256": sha256(path),
                          "n_proteins": df["name"].nunique(), "contrasts": ", ".join(contrasts)}])
    return LoadedInputs(long, contrasts, prov)
