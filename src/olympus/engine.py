"""Per-temperature moderated linear models, with Deimos' limma engine.

Olympus does not carry its own copy of limma: the empirical Bayes engine is
Deimos' ``limma_ebayes.py`` (validated against R limma 3.58), imported from the
Deimos folder. This keeps one statistical engine for the platform.

Where to find it, in order: ``deimos_dir`` argument (config ``shift_test.deimos_dir``),
environment variable ``DEIMOS_DIR``, or a ``limma_ebayes`` module already importable.

Paired design
-------------
In a TPP experiment one lysate (one biological replicate) is usually split
into every condition and every temperature. The replicate is then a block:
the per-temperature model becomes ``~0 + condition + block`` (Deimos' paired
design), which removes the lysate effect from the residual variance. With
``paired="auto"`` the block is used when every block label appears in every
condition.
"""

from __future__ import annotations

import importlib
import logging
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

_LIMMA = None


def load_limma(deimos_dir: str | Path | None = None):
    """Import Deimos' ``limma_ebayes`` module (cached)."""
    global _LIMMA
    if _LIMMA is not None and deimos_dir is None:
        return _LIMMA
    tried = []
    for cand in (deimos_dir, os.environ.get("DEIMOS_DIR")):
        if not cand:
            continue
        p = Path(cand).expanduser().resolve()
        if p.is_file():
            p = p.parent
        tried.append(str(p))
        if (p / "limma_ebayes.py").exists():
            if str(p) not in sys.path:
                sys.path.insert(0, str(p))
            break
    try:
        mod = importlib.import_module("limma_ebayes")
    except ImportError as e:
        where = f" (looked in: {', '.join(tried)})" if tried else ""
        raise ImportError(
            "The shift test needs Deimos' limma engine (limma_ebayes.py). Set shift_test.deimos_dir "
            "in the YAML, or the DEIMOS_DIR environment variable, to your Deimos folder" + where) from e
    for fn in ("lm_fit", "contrasts_fit", "ebayes"):
        if not hasattr(mod, fn):
            raise ImportError(f"{getattr(mod, '__file__', mod)} has no {fn}(): not Deimos' limma_ebayes")
    _LIMMA = mod
    return mod


def parse_contrast(name: str) -> tuple[str, str]:
    if "_vs_" not in name:
        raise ValueError(f"Contrast '{name}' is not of the form A_vs_B")
    a, b = name.split("_vs_", 1)
    return a, b


@dataclass
class TemperatureFit:
    temperature: str
    proteins: pd.Index
    diff: np.ndarray          # n × n_contrasts (log2 FC)
    se: np.ndarray            # n × n_contrasts (moderated standard error)
    p: np.ndarray             # n × n_contrasts
    df_total: float
    resid: np.ndarray         # n × n_samples (after condition [+ block] effects)
    units: list[str]          # 'condition|block' per residual column
    paired: bool
    d0: float
    contrasts: list[str] = field(default_factory=list)


def complete_blocks(conditions: list[str], blocks: list[str]) -> bool:
    groups = list(dict.fromkeys(conditions))
    by_block = {}
    for c, b in zip(conditions, blocks):
        by_block.setdefault(b, set()).add(c)
    return len(by_block) >= 2 and all(set(groups) <= v for v in by_block.values())


def design_matrix(conditions: list[str], blocks: list[str] | None,
                  groups: list[str]) -> np.ndarray:
    cols = [[1.0 if c == g else 0.0 for g in groups] for c in conditions]
    D = np.array(cols, float)
    if blocks is not None:
        levels = list(dict.fromkeys(blocks))[1:]        # first block = reference
        B = np.array([[1.0 if b == lv else 0.0 for lv in levels] for b in blocks], float)
        D = np.hstack([D, B])
    return D


def fit_temperature(limma, temperature: str, mat: pd.DataFrame, conditions: list[str],
                    blocks: list[str], contrasts: list[str], paired: str | bool = "auto") -> TemperatureFit:
    """limma fit of one temperature: ~0 + condition (+ block when paired)."""
    use_block = False
    if paired in (True, "true", "True") or paired == "auto":
        ok = complete_blocks(conditions, blocks)
        if paired == "auto":
            use_block = ok
        elif not ok:
            raise ValueError(f"[{temperature}] paired fit requested but the blocks are incomplete")
        else:
            use_block = True
    groups = list(dict.fromkeys(conditions))
    D = design_matrix(conditions, blocks if use_block else None, groups)
    C = np.zeros((D.shape[1], len(contrasts)))
    for j, name in enumerate(contrasts):
        a, b = parse_contrast(name)
        if a not in groups or b not in groups:
            raise ValueError(f"[{temperature}] contrast {name}: conditions {groups}")
        C[groups.index(a), j] = 1.0
        C[groups.index(b), j] = -1.0
    X = mat.to_numpy(float)
    fit = limma.lm_fit(X, D)
    eb = limma.ebayes(limma.contrasts_fit(fit, C))
    se = np.asarray(eb["stdev_unscaled"]) * np.sqrt(np.asarray(eb["s2_post"]))[:, None]
    resid = X - np.asarray(fit["coefficients"]) @ D.T
    units = [f"{c}|{b}" for c, b in zip(conditions, blocks)]
    return TemperatureFit(temperature, mat.index, np.asarray(eb["coefficients"]), se,
                          np.asarray(eb["p_value"]), float(eb["df_total"]), resid, units,
                          use_block, float(eb.get("d0", np.nan)), list(contrasts))


def fit_series(limma, matrices: dict, temperatures: list[str], contrasts: list[str],
               paired: str | bool = "auto", relabel: dict | None = None,
               proteins: set | None = None) -> dict[str, TemperatureFit]:
    """Fit every temperature. ``relabel`` maps a unit 'condition|block' to a new
    condition (used by the permutations); ``proteins`` restricts the universe."""
    fits = {}
    for t in temperatures:
        tm = matrices[t]
        mat = tm.imputed
        if proteins is not None:
            mat = mat[mat.index.isin(proteins)]
        conds = tm.samples["condition"].astype(str).tolist()
        blocks = tm.block_labels()
        if relabel:
            conds = [relabel.get(f"{c}|{b}", c) for c, b in zip(conds, blocks)]
        fits[t] = fit_temperature(limma, t, mat, conds, blocks, contrasts, paired)
    return fits
