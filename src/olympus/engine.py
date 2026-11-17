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
design), which removes the lysate effect from the residual variance. The block
costs degrees of freedom (n_blocks − 1 per temperature), so it only pays when the
lysate effect is real. With ``paired="auto"`` the block is used when (1) every
block label appears in every condition at every temperature and (2) the data show
a positive within-lysate correlation (``block_correlation``).

Block correlation
-----------------
At each temperature, residuals from the condition means are standardised per
protein; for every pair of samples from different conditions, the mean product
over proteins estimates their correlation. Same-lysate pairs estimate the
intra-lysate correlation ρ; different-lysate pairs estimate −ρ/(n−1) (n blocks),
and both are 0 without a lysate effect. δ_T = mean(same) − mean(different) is
tested > 0 over the temperatures (one-sided t-test). ``paired="auto"`` pairs when
p < 0.05.
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


@dataclass
class BlockCorrelation:
    rho: float                    # mean same-lysate correlation (estimate of ρ)
    different: float              # mean different-lysate correlation (≈ −ρ/(n−1))
    delta: pd.Series              # per temperature: same − different
    t: float
    p: float
    complete: bool                # blocks complete at every temperature
    n_blocks: int

    def decide(self, alpha: float = 0.05) -> bool:
        return bool(self.complete and np.isfinite(self.p) and self.p < alpha and self.rho > 0)

    def describe(self) -> str:
        if not self.complete:
            return "blocks incomplete (some lysate missing a condition)"
        return (f"within-lysate correlation {self.rho:.2f} (different lysates {self.different:+.2f}); "
                f"δ > 0 at {int((self.delta > 0).sum())}/{len(self.delta)} temperatures, "
                f"one-sided t = {self.t:.2f}, p = {self.p:.2g}")


def block_correlation(matrices: dict, temperatures: list[str],
                      proteins: set | None = None) -> BlockCorrelation:
    """Within-lysate residual correlation between conditions (see module docstring)."""
    from scipy import stats

    deltas, same_all, diff_all, complete = {}, [], [], True
    n_blocks = 0
    for t in temperatures:
        tm = matrices[t]
        mat = tm.imputed if proteins is None else tm.imputed[tm.imputed.index.isin(proteins)]
        conds = tm.samples["condition"].astype(str).tolist()
        blocks = tm.block_labels()
        if not complete_blocks(conds, blocks):
            complete = False
            continue
        n_blocks = max(n_blocks, len(set(blocks)))
        X = mat.to_numpy(float)
        E = X.copy()
        for g in set(conds):
            idx = [i for i, c in enumerate(conds) if c == g]
            E[:, idx] -= X[:, idx].mean(1, keepdims=True)
        rms = np.sqrt((E ** 2).mean(1))
        Z = E[rms > 0] / rms[rms > 0, None]
        same, diff = [], []
        for i in range(len(conds)):
            for j in range(i + 1, len(conds)):
                if conds[i] == conds[j]:
                    continue
                c = float(np.mean(Z[:, i] * Z[:, j]))
                (same if blocks[i] == blocks[j] else diff).append(c)
        if same and diff:
            deltas[t] = float(np.mean(same) - np.mean(diff))
            same_all.append(np.mean(same))
            diff_all.append(np.mean(diff))
    d = pd.Series(deltas, dtype=float)
    if len(d) >= 3 and d.std() > 0:
        res = stats.ttest_1samp(d.to_numpy(), 0.0, alternative="greater")
        tval, pval = float(res.statistic), float(res.pvalue)
    else:
        tval, pval = np.nan, np.nan
    return BlockCorrelation(float(np.mean(same_all)) if same_all else np.nan,
                            float(np.mean(diff_all)) if diff_all else np.nan,
                            d, tval, pval, complete and len(d) > 0, n_blocks)


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
