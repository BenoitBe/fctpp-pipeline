"""Shift test: protein-level, FDR-controlled detection of temperature-localised
condition effects, plus permutation estimates of false discoveries.

Why
---
The legacy rule ("|log2FC| ≥ 0.6 and p < 0.01 at 2 adjacent temperatures")
thresholds raw per-temperature p-values: about 2,000 proteins × 12 temperatures
are screened per contrast with no multiplicity control, and adjacent
temperatures share the same lysates. On the TPP data, within-lysate permutations
give about as many "hits" as the real labels.

Statistic
---------
Per temperature T and contrast A_vs_B, Deimos' limma (paired on the lysate when
possible) gives Δ_T and its moderated standard error SE_T. For each window of
``w`` adjacent temperatures:

    D_w = mean(Δ_T, T ∈ w) − baseline

* ``baseline: median`` (default): median of the protein's Δ profile over all its
  temperatures. A thermal shift moves a few temperatures, so the median stays at
  the protein's baseline; a constant offset (handling artefact, or abundance
  change in cells) is removed. Var ≈ (π/2)·mean(SE²)/K.
* ``baseline: lowest``: Δ at the lowest temperature of the series (classical
  ΔΔ / interaction contrast; noisier: one estimate).
* ``baseline: none``: Δ itself (lysate experiments with no offset expected).

Var(D_w) = [Σ SE_T² + 2 Σ ρ_{T,T+1} SE_T SE_{T+1}] / w² + Var(baseline), where
ρ is the consensus correlation of residuals between adjacent temperatures for
the same lysate (estimated across proteins, floored at 0). Welch–Satterthwaite
degrees of freedom; two-sided t p-value per window.

Protein level: the window p-values are combined with Simes (valid under the
positive dependence of overlapping windows), then Benjamini–Hochberg across
proteins, per contrast. A protein is a shift hit when q ≤ alpha and the best
window has |D_w| ≥ fc_cutoff.

Permutations
------------
Under H0 (no A–B difference), A and B are exchangeable *within a lysate*. For a
contrast, the labels A/B are swapped in a subset of lysates, the same subset at
every temperature (keeping the between-temperature correlation), excluding the
identity and the full swap (which only flips signs). Each pattern gives a full
null dataset; the mean number of hits estimates the expected false discoveries
of a rule (FDR ≈ mean null hits / observed hits). With 3 lysates there are only
6 such patterns: the estimate is indicative. Partial swaps keep part of a real
signal, so it tends to be conservative.
"""

from __future__ import annotations

import itertools
import logging
import warnings
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
from scipy import stats

from .engine import TemperatureFit, fit_series, parse_contrast

log = logging.getLogger(__name__)

BASELINES = ("median", "lowest", "none")


@dataclass
class ShiftParams:
    baseline: str = "median"
    window: int = 2
    alpha: float = 0.05
    fc_cutoff: float = 0.6
    use_adjacent_correlation: bool = True
    paired: str = "auto"

    def validate(self) -> None:
        if self.baseline not in BASELINES:
            raise ValueError(f"shift_test.baseline must be one of {BASELINES}")
        if self.window < 1:
            raise ValueError("shift_test.window must be >= 1")
        if not 0 < self.alpha < 1:
            raise ValueError("shift_test.alpha must be in (0, 1)")


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg (p.adjust 'BH'); NaN stay NaN and are not counted."""
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    n = int(ok.sum())
    if n == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(n)
    res[order] = np.minimum(ranked, 1.0)
    out[ok] = res
    return out


def simes(p_rows: np.ndarray) -> np.ndarray:
    """Simes combination per row (NaN ignored; all-NaN rows give NaN)."""
    out = np.full(p_rows.shape[0], np.nan)
    for i, row in enumerate(p_rows):
        pv = np.sort(row[np.isfinite(row)])
        m = len(pv)
        if m:
            out[i] = min(1.0, float(np.min(pv * m / np.arange(1, m + 1))))
    return out


def adjacent_correlation(fits: dict[str, TemperatureFit], temperatures: list[str],
                         floor: float = 0.0, cap: float = 0.95) -> dict[tuple[str, str], float]:
    """Consensus correlation (Fisher-z mean over proteins) of the residuals of the
    same unit (condition|lysate) at adjacent temperatures."""
    out = {}
    for a, b in zip(temperatures[:-1], temperatures[1:]):
        fa, fb = fits[a], fits[b]
        common_u = [u for u in fa.units if u in fb.units]
        common_p = fa.proteins.intersection(fb.proteins)
        if len(common_u) < 3 or len(common_p) < 10:
            out[(a, b)] = 0.0
            continue
        ra = fa.resid[fa.proteins.get_indexer(common_p)][:, [fa.units.index(u) for u in common_u]]
        rb = fb.resid[fb.proteins.get_indexer(common_p)][:, [fb.units.index(u) for u in common_u]]
        num = (ra * rb).sum(1)
        den = np.sqrt((ra ** 2).sum(1) * (rb ** 2).sum(1))
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.clip(num / den, -0.999, 0.999)
        r = r[np.isfinite(r)]
        rho = float(np.tanh(np.mean(np.arctanh(r)))) if len(r) else 0.0
        out[(a, b)] = float(min(max(rho, floor), cap))
    return out


def _aligned(fits, temperatures, names, attr, j):
    cols = []
    for t in temperatures:
        f = fits[t]
        s = pd.Series(getattr(f, attr)[:, j], index=f.proteins)
        cols.append(s.reindex(names).to_numpy(float))
    return np.column_stack(cols)


def shift_test(fits: dict[str, TemperatureFit], temperatures: list[str], contrasts: list[str],
               params: ShiftParams, rho: dict | None = None,
               universe: list[str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Protein-level table (one row per contrast × protein) and window table."""
    params.validate()
    names = universe or sorted(set().union(*[set(fits[t].proteins) for t in temperatures]))
    names = list(names)
    K, w = len(temperatures), params.window
    dfs = np.array([fits[t].df_total for t in temperatures])
    rho = rho or {}
    adj = np.array([rho.get((a, b), 0.0) for a, b in zip(temperatures[:-1], temperatures[1:])]) \
        if params.use_adjacent_correlation else np.zeros(max(K - 1, 0))
    prot_rows, win_rows = [], []
    for j, c in enumerate(contrasts):
        D = _aligned(fits, temperatures, names, "diff", j)
        SE = _aligned(fits, temperatures, names, "se", j)
        present = np.isfinite(D)
        kp = present.sum(1)
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)     # all-NaN rows -> NaN
            if params.baseline == "median":
                base = np.nanmedian(D, axis=1)
                var_b = (np.pi / 2) * np.nanmean(SE ** 2, axis=1) / np.maximum(kp, 1)
                df_b = np.maximum(kp, 1) * float(np.mean(dfs))
            elif params.baseline == "lowest":
                base, var_b = D[:, 0], SE[:, 0] ** 2
                df_b = np.full(len(names), dfs[0])
            else:
                base, var_b = np.zeros(len(names)), np.zeros(len(names))
                df_b = np.full(len(names), np.inf)
        starts = [i for i in range(K - w + 1) if not (params.baseline == "lowest" and i == 0)]
        P = np.full((len(names), len(starts)), np.nan)
        DW = np.full_like(P, np.nan)
        SW = np.full_like(P, np.nan)
        for k, i in enumerate(starts):
            idx = list(range(i, i + w))
            dw = D[:, idx].mean(1) - base                       # NaN if any temperature is missing
            v = (SE[:, idx] ** 2).sum(1)
            for a in range(i, i + w - 1):
                v = v + 2 * adj[a] * SE[:, a] * SE[:, a + 1]
            vw = v / w ** 2 + var_b
            comps = [((SE[:, a] ** 2) / w ** 2, dfs[a]) for a in idx]
            num = sum(x for x, _ in comps) + var_b
            den = sum(x ** 2 / d for x, d in comps) + np.where(np.isfinite(df_b), var_b ** 2 / df_b, 0.0)
            with np.errstate(invalid="ignore", divide="ignore"):
                dfw = num ** 2 / den
                tw = dw / np.sqrt(vw)
                pw = 2 * stats.t.sf(np.abs(tw), dfw)
            P[:, k], DW[:, k], SW[:, k] = pw, dw, np.sqrt(vw)
            lab = f"{temperatures[i]}-{temperatures[i + w - 1]}" if w > 1 else temperatures[i]
            ok = np.isfinite(pw)
            win_rows.append(pd.DataFrame({"contrast": c, "name": np.array(names)[ok], "window": lab,
                                          "dd": dw[ok], "se": np.sqrt(vw)[ok], "df": dfw[ok], "p": pw[ok]}))
        p_prot = simes(P)
        q = bh(p_prot)
        anyp = np.isfinite(P).any(1)
        best = np.where(anyp, np.nanargmin(np.where(np.isfinite(P), P, np.inf), axis=1), 0)
        rows = np.arange(len(names))
        dd_best = np.where(anyp, DW[rows, best], np.nan)
        lab_best = [(f"{temperatures[starts[b]]}-{temperatures[starts[b] + w - 1]}" if w > 1
                     else temperatures[starts[b]]) if anyp[i] else None for i, b in enumerate(best)]
        tab = pd.DataFrame({
            "contrast": c, "name": names, "n_temperatures": kp, "baseline": base,
            "best_window": lab_best, "dd_best": dd_best, "se_best": np.where(anyp, SW[rows, best], np.nan),
            "p_best_window": np.where(anyp, P[rows, best], np.nan), "n_windows": np.isfinite(P).sum(1),
            "p_protein": p_prot, "q_protein": q,
        })
        tab["direction"] = np.where(tab["dd_best"] > 0, "up", np.where(tab["dd_best"] < 0, "down", None))
        tab["hit"] = (tab["q_protein"] <= params.alpha) & (tab["dd_best"].abs() >= params.fc_cutoff)
        prot_rows.append(tab)
    prot = pd.concat(prot_rows, ignore_index=True)
    wins = pd.concat(win_rows, ignore_index=True) if win_rows else pd.DataFrame()
    return prot, wins


# ----------------------------------------------------------------------------- permutations

def swap_patterns(blocks: list[str], max_patterns: int = 64, seed: int = 42) -> list[tuple[int, ...]]:
    """Subsets of lysates in which A and B are swapped (identity and full swap excluded)."""
    pats = [p for p in itertools.product([0, 1], repeat=len(blocks)) if 0 < sum(p) < len(blocks)]
    # a pattern and its complement give the same |statistics| (sign flip): keep one of each
    seen, uniq = set(), []
    for p in pats:
        comp = tuple(1 - x for x in p)
        if comp in seen:
            continue
        seen.add(p)
        uniq.append(p)
    if len(uniq) > max_patterns:
        rng = np.random.default_rng(seed)
        uniq = [uniq[i] for i in sorted(rng.choice(len(uniq), size=max_patterns, replace=False))]
    return uniq


def permutation_null(limma, matrices: dict, temperatures: list[str], contrasts: list[str],
                     evaluate, paired="auto", max_patterns: int = 64, seed: int = 42,
                     proteins: set | None = None) -> dict[str, list[int]]:
    """Null hit counts per contrast. ``evaluate(fits, contrast) -> int``."""
    blocks = sorted(set(matrices[temperatures[0]].block_labels()))
    out = {}
    for c in contrasts:
        a, b = parse_contrast(c)
        counts = []
        for pat in swap_patterns(blocks, max_patterns, seed):
            relabel = {}
            for blk, flip in zip(blocks, pat):
                if flip:
                    relabel[f"{a}|{blk}"] = b
                    relabel[f"{b}|{blk}"] = a
            fits = fit_series(limma, matrices, temperatures, [c], paired, relabel, proteins)
            counts.append(int(evaluate(fits, c)))
        out[c] = counts
        log.info("Permutations %s: %d patterns, null hits %s", c, len(counts), counts)
    return out


def fits_to_long(fits: dict[str, TemperatureFit], temperatures: list[str], contrasts: list[str]) -> pd.DataFrame:
    """Long table (Temperature, name, <c>_diff, <c>_p.val) from refits, for the legacy rule."""
    parts = []
    for t in temperatures:
        f = fits[t]
        d = {"Temperature": t, "name": f.proteins.astype(str)}
        for j, c in enumerate(f.contrasts):
            if c in contrasts:
                d[f"{c}_diff"] = f.diff[:, j]
                d[f"{c}_p.val"] = f.p[:, j]
                d[f"{c}_p.adj"] = bh(f.p[:, j])
        parts.append(pd.DataFrame(d))
    return pd.concat(parts, ignore_index=True)


def params_dict(p: ShiftParams) -> dict:
    return asdict(p)
