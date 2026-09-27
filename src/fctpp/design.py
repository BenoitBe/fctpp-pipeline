"""Experimental design reading and sample ↔ column matching.

R behaviour being replaced
--------------------------
``find_design()`` read one Excel design whose labels are those of the *first*
temperature (``Sple.1B``, ``Sple.1C``...). ``DEP::make_se()`` then matched it
against every temperature (``Sple.2B``..., ``Sple.12J``) because it strips the
longest common prefix on both sides (``delete_prefix``) before matching, so only
the trailing letter (B..J) was actually compared. That worked by side effect and
was never reported.

Here the two strategies are explicit, and the resulting mapping is written to
``sample_mapping.csv`` for every temperature:

* ``exact``  — labels equal the cleaned column names (after canonicalisation).
* ``suffix`` — DEP-compatible: match on what remains after removing the longest
  common prefix of each side; the match must be a bijection.

A *long* design (with a ``temperature`` column) gives explicit labels per
temperature and is the recommended format.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pandas as pd

# Deimos cleans DIA-NN run names with this regex (deimos.py, load_data);
# the per-temperature design must use the same cleaned labels.
DEIMOS_SAMPLE_REGEX = r"(Sple-[^_]+)_"

REQUIRED = ("label", "condition", "replicate")


def canonical(label: object) -> str:
    """R ``make.names``-like canonical form: every character outside
    [A-Za-z0-9._] becomes '.', so 'Sple-1B' == 'Sple.1B' (R check.names)."""
    return re.sub(r"[^A-Za-z0-9._]", ".", str(label).strip())


def clean_sample_name(column: str, regex: str | None = DEIMOS_SAMPLE_REGEX) -> str:
    """Sample label derived from a raw DIA-NN run column, exactly as Deimos does
    (first capture group of ``regex``; the full name if it does not match)."""
    if not regex:
        return str(column)
    m = re.search(regex, str(column))
    return m.group(1) if m else str(column)


def longest_common_prefix(values: list[str]) -> str:
    if not values:
        return ""
    return os.path.commonprefix([str(v) for v in values])


def read_design(path: str | Path, sheet: int | str = 0) -> pd.DataFrame:
    """Read a design table (.xlsx/.xls/.csv/.tsv; ';', ',' or tab separated)."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Design file not found: {path}")
    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path, sheet_name=sheet, dtype=str)
    else:
        df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    df.columns = [str(c).strip() for c in df.columns]
    lower = {c.lower(): c for c in df.columns}
    df = df.rename(columns={lower[k]: k for k in list(REQUIRED) + ["temperature"] if k in lower})
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing design column(s) {missing}; found {list(df.columns)}")
    for c in df.columns:
        df[c] = df[c].astype("string").str.strip()
    df = df[df["label"].notna() & (df["label"] != "")].reset_index(drop=True)
    df["replicate"] = df["replicate"].str.replace(r"\.0$", "", regex=True)
    key = ["temperature", "label"] if "temperature" in df.columns else ["label"]
    dup = df.duplicated(key)
    if dup.any():
        raise ValueError(f"{path}: duplicated label(s): {df.loc[dup, 'label'].tolist()}")
    return df.astype(object).where(df.notna(), None)


def find_design(folder: str | Path, design_name: str = "experimental_design") -> Path:
    """R ``find_design()``: first file in ``folder`` whose name starts with
    ``design_name`` (case-insensitive — R's default folder 'tpp_data' did not
    exist on case-sensitive file systems where the data sat in 'TPP_data')."""
    folder = Path(folder)
    if not folder.is_dir():
        # case-insensitive fallback on the folder name itself
        parent = folder.parent if str(folder.parent) else Path(".")
        cands = [p for p in parent.iterdir() if p.is_dir() and p.name.lower() == folder.name.lower()] \
            if parent.is_dir() else []
        if not cands:
            raise FileNotFoundError(f"Folder not found: {folder}")
        folder = cands[0]
    pat = re.compile(re.escape(design_name) + r".*\.(xlsx?|csv|tsv)$", re.IGNORECASE)
    files = sorted(p for p in folder.iterdir() if pat.match(p.name))
    if not files:
        raise FileNotFoundError(f"No design file '{design_name}*' in {folder}")
    return files[0]


def match_samples(columns: list[str], design: pd.DataFrame, *,
                  sample_regex: str | None = DEIMOS_SAMPLE_REGEX,
                  strategy: str = "auto") -> pd.DataFrame:
    """Map design rows onto raw sample columns of one temperature.

    Returns one row per design sample: ``column`` (raw name), ``sample_label``
    (cleaned name, used as label in the per-temperature Deimos design),
    ``design_label``, ``condition``, ``replicate``, ``match``.
    """
    cleaned = [clean_sample_name(c, sample_regex) for c in columns]
    if len(set(cleaned)) != len(cleaned):
        dups = sorted({x for x in cleaned if cleaned.count(x) > 1})
        raise ValueError(f"Several columns give the same sample label after cleaning: {dups}")
    canon_cols = [canonical(c) for c in cleaned]
    labels = [str(x) for x in design["label"]]
    canon_lab = [canonical(x) for x in labels]

    def _result(idx: list[int], how: str) -> pd.DataFrame:
        return pd.DataFrame({
            "design_label": labels,
            "column": [columns[i] for i in idx],
            "sample_label": [cleaned[i] for i in idx],
            "condition": list(design["condition"]),
            "replicate": list(design["replicate"]),
            "match": how,
        })

    if strategy in ("auto", "exact"):
        lookup = {c: i for i, c in enumerate(canon_cols)}
        idx = [lookup.get(x) for x in canon_lab]
        if all(i is not None for i in idx):
            return _result(idx, "exact")
        if strategy == "exact":
            miss = [lab for lab, i in zip(labels, idx) if i is None]
            raise ValueError(f"Design labels without an exact column match: {miss}; "
                             f"columns (cleaned): {cleaned}")

    if strategy in ("auto", "suffix"):
        p_lab = longest_common_prefix(canon_lab)
        p_col = longest_common_prefix(canon_cols)
        suf_lab = [x[len(p_lab):] for x in canon_lab]
        suf_col = [x[len(p_col):] for x in canon_cols]
        if len(set(suf_lab)) != len(suf_lab) or "" in suf_lab:
            raise ValueError(f"Suffix matching impossible: design suffixes not unique ({suf_lab})")
        lookup = {}
        for i, s in enumerate(suf_col):
            lookup.setdefault(s, []).append(i)
        idx, problems = [], []
        for lab, s in zip(labels, suf_lab):
            hit = lookup.get(s, [])
            if len(hit) != 1:
                problems.append(f"{lab} -> suffix '{s}' matches {len(hit)} column(s)")
                idx.append(None)
            else:
                idx.append(hit[0])
        if problems:
            raise ValueError("Suffix matching failed (prefixes '%s' / '%s'): %s"
                             % (p_lab, p_col, "; ".join(problems)))
        return _result(idx, f"suffix (prefix design='{p_lab}', columns='{p_col}')")

    raise ValueError(f"Unknown matching strategy {strategy!r}")
