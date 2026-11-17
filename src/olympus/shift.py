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

* ``baseline: none`` (default): D_w = mean Δ over the window. Most powerful;
  combined with the localisation check below it does not call constant offsets.
* ``baseline: median``: median of the protein's Δ profile over all its
  temperatures. A thermal shift moves a few temperatures, so the median stays at
  the protein's baseline; a constant offset (handling artefact, or abundance
  change in cells) is removed. Var ≈ (π/2)·mean(SE²)/K. Less powerful.
* ``baseline: lowest``: Δ at the lowest temperature of the series (classical
  ΔΔ / interaction contrast; noisiest: one estimate).

Localisation check (``localisation_p``, default 0.01, with ``none``/``lowest``):
the detected window must also differ from the median of the profile (window
p-value ≤ localisation_p). A constant offset is significant everywhere but never
localised, so it is not called a shift.

Validation (docs/validation/simulate_shift.py, realistic noise: protein-specific
variances d0 = 2, t5 noise, censored + imputed low values; 3 lysates, windows
2/3/4): false discovery proportion ≈ 2 %, sensitivity 0.75 for ΔTm = 2 °C and
≥ 0.93 for ΔTm ≥ 4 °C, < 1 of 40 constant offsets called (the legacy rule calls
~33 of 40).

Var(D_w) = [Σ SE_T² + 2 Σ ρ_{T,T+1} SE_T SE_{T+1}] / w² + Var(baseline), where
ρ is the consensus correlation of residuals between adjacent temperatures for
the same lysate (estimated across proteins, floored at 0). Welch–Satterthwaite
degrees of freedom; two-sided t p-value per window.

Several window sizes (``windows: [2, 3, 4]``): a narrow shift is best seen on 2
temperatures, a broad one (large ΔTm, shallow melting curve) on 3–4; every window
of every size enters the same Simes combination.

Protein level: the window p-values are combined with Simes (valid under the
positive dependence of overlapping windows), then Benjamini–Hochberg across
proteins, per contrast. A protein is a shift hit when q ≤ alpha and the effect
size ``dd_peak`` reaches fc_cutoff: the largest same-sign |D| among the windows
lying inside the best window (= |D| of the best window with a single size).

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
    baseline: str = "none"
    window: int | None = None                # single window size; when set, overrides ``windows``
    windows: tuple[int, ...] = (2, 3, 4)     # window sizes combined by Simes
    alpha: float = 0.05
    fc_cutoff: float = 0.6
    use_adjacent_correlation: bool = True
    paired: str = "auto"
    localisation_p: float | None = 0.01   # with baseline none/lowest; None disables the check

    @property
    def window_sizes(self) -> tuple[int, ...]:
        if self.window is not None:
            return (int(self.window),)
        return tuple(sorted({int(w) for w in self.windows})) if self.windows else (2,)

    def describe_windows(self) -> str:
        ws = self.window_sizes
        return f"window {ws[0]}" if len(ws) == 1 else "windows " + "/".join(map(str, ws))

    def validate(self) -> None:
        if self.baseline not in BASELINES:
            raise ValueError(f"shift_test.baseline must be one of {BASELINES}")
        if (self.window is not None and int(self.window) < 1) or any(int(w) < 1 for w in (self.windows or ())):
            raise ValueError("shift_test.window(s) must be >= 1")
        if not 0 < self.alpha < 1:
            raise ValueError("shift_test.alpha must be in (0, 1)")
        if self.localisation_p is not None and not 0 < self.localisation_p <= 1:
            raise ValueError("shift_test.localisation_p must be in (0, 1] or null")


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


def _baseline(D, SE, dfs, kind):
    n = D.shape[0]
    kp = np.isfinite(D).sum(1)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)          # all-NaN rows -> NaN
        if kind == "median":
            return (np.nanmedian(D, axis=1),
                    (np.pi / 2) * np.nanmean(SE ** 2, axis=1) / np.maximum(kp, 1),
                    np.maximum(kp, 1) * float(np.mean(dfs)))
        if kind == "lowest":
            return D[:, 0], SE[:, 0] ** 2, np.full(n, dfs[0])
    return np.zeros(n), np.zeros(n), np.full(n, np.inf)


def _window_stats(D, SE, dfs, adj, starts, w, base, var_b, df_b):
    """D, SE, df and p of every window (columns), Welch–Satterthwaite df."""
    shape = (D.shape[0], len(starts))
    P, DW, SW, DF = (np.full(shape, np.nan) for _ in range(4))
    for k, i in enumerate(starts):
        idx = list(range(i, i + w))
        dw = D[:, idx].mean(1) - base                           # NaN if a temperature is missing
        v = (SE[:, idx] ** 2).sum(1)
        for a in range(i, i + w - 1):
            v = v + 2 * adj[a] * SE[:, a] * SE[:, a + 1]
        vw = v / w ** 2 + var_b
        comps = [((SE[:, a] ** 2) / w ** 2, dfs[a]) for a in idx]
        num = sum(x for x, _ in comps) + var_b
        den = sum(x ** 2 / d for x, d in comps) + np.where(np.isfinite(df_b), var_b ** 2 / df_b, 0.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            dfw = num ** 2 / den
            pw = 2 * stats.t.sf(np.abs(dw / np.sqrt(vw)), dfw)
        P[:, k], DW[:, k], SW[:, k], DF[:, k] = pw, dw, np.sqrt(vw), dfw
    return P, DW, SW, DF


def _window_grid(temperatures, sizes, baseline):
    """(size, start, label) of every window, sizes in increasing order."""
    K = len(temperatures)
    grid = []
    for w in sizes:
        for i in range(K - w + 1):
            if baseline == "lowest" and i == 0:
                continue
            lab = f"{temperatures[i]}-{temperatures[i + w - 1]}" if w > 1 else temperatures[i]
            grid.append((w, i, lab))
    return grid


def _all_windows(D, SE, dfs, adj, grid, base, var_b, df_b):
    """Window statistics for a grid mixing several sizes (columns in grid order)."""
    shape = (D.shape[0], len(grid))
    P, DW, SW, DF = (np.full(shape, np.nan) for _ in range(4))
    for w in sorted({g[0] for g in grid}):
        cols = [k for k, g in enumerate(grid) if g[0] == w]
        starts = [grid[k][1] for k in cols]
        p, dw, sw, df = _window_stats(D, SE, dfs, adj, starts, w, base, var_b, df_b)
        P[:, cols], DW[:, cols], SW[:, cols], DF[:, cols] = p, dw, sw, df
    return P, DW, SW, DF


def _peak_within(DW, grid, best, sign):
    """Largest same-sign |D| among the windows (any size) lying inside each protein's
    best window, the best window included: the effect size the fc cutoff is applied
    to (= D of the best window when a single size is used)."""
    out = np.full(DW.shape[0], np.nan)
    spans = [(g[1], g[1] + g[0]) for g in grid]
    for i, b in enumerate(best):
        if not np.isfinite(sign[i]) or sign[i] == 0:
            continue
        sb, eb = spans[b]
        inside = [k for k, (s0, e0) in enumerate(spans) if s0 >= sb and e0 <= eb]
        vals = DW[i, inside] * sign[i]
        vals = vals[np.isfinite(vals)]
        out[i] = sign[i] * vals.max() if len(vals) else DW[i, b]
    return out


def shift_test(fits: dict[str, TemperatureFit], temperatures: list[str], contrasts: list[str],
               params: ShiftParams, rho: dict | None = None,
               universe: list[str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Protein-level table (one row per contrast × protein) and window table."""
    params.validate()
    names = list(universe or sorted(set().union(*[set(fits[t].proteins) for t in temperatures])))
    K = len(temperatures)
    sizes = tuple(w for w in params.window_sizes if w <= K - (1 if params.baseline == "lowest" else 0))
    if not sizes:
        raise ValueError(f"shift_test: no window size fits {K} temperatures")
    dfs = np.array([fits[t].df_total for t in temperatures])
    rho = rho or {}
    adj = np.array([rho.get((a, b), 0.0) for a, b in zip(temperatures[:-1], temperatures[1:])]) \
        if params.use_adjacent_correlation else np.zeros(max(K - 1, 0))
    grid = _window_grid(temperatures, sizes, params.baseline)
    labels = [g[2] for g in grid]
    check_loc = params.localisation_p is not None and params.baseline != "median"
    prot_rows, win_rows = [], []
    rows = np.arange(len(names))
    for j, c in enumerate(contrasts):
        D = _aligned(fits, temperatures, names, "diff", j)
        SE = _aligned(fits, temperatures, names, "se", j)
        kp = np.isfinite(D).sum(1)
        base, var_b, df_b = _baseline(D, SE, dfs, params.baseline)
        P, DW, SW, DF = _all_windows(D, SE, dfs, adj, grid, base, var_b, df_b)
        for k, (w, _, lab) in enumerate(grid):
            ok = np.isfinite(P[:, k])
            win_rows.append(pd.DataFrame({"contrast": c, "name": np.array(names)[ok], "window": lab,
                                          "size": w, "dd": DW[ok, k], "se": SW[ok, k], "df": DF[ok, k],
                                          "p": P[ok, k]}))
        p_prot = simes(P)
        q = bh(p_prot)
        anyp = np.isfinite(P).any(1)
        best = np.where(anyp, np.argmin(np.where(np.isfinite(P), P, np.inf), axis=1), 0)
        dd_best = np.where(anyp, DW[rows, best], np.nan)
        dd_peak = _peak_within(DW, grid, best, np.sign(dd_best))
        tab = pd.DataFrame({
            "contrast": c, "name": names, "n_temperatures": kp, "baseline": base,
            "best_window": [labels[b] if anyp[i] else None for i, b in enumerate(best)],
            "window_size": np.where(anyp, [grid[b][0] for b in best], 0),
            "dd_best": dd_best, "dd_peak": np.where(anyp, dd_peak, np.nan),
            "se_best": np.where(anyp, SW[rows, best], np.nan),
            "p_best_window": np.where(anyp, P[rows, best], np.nan), "n_windows": np.isfinite(P).sum(1),
            "p_protein": p_prot, "q_protein": q,
        })
        tab["direction"] = np.where(tab["dd_best"] > 0, "up", np.where(tab["dd_best"] < 0, "down", None))
        hit = (tab["q_protein"] <= params.alpha) & (tab["dd_peak"].abs() >= params.fc_cutoff)
        if check_loc:
            # localisation: the detected window must also stand out from the profile median,
            # which rejects constant offsets (handling artefacts, abundance changes)
            mb, mv, mdf = _baseline(D, SE, dfs, "median")
            Pm, DWm, _, _ = _all_windows(D, SE, dfs, adj, grid, mb, mv, mdf)
            tab["median_baseline"] = mb
            tab["dd_vs_median"] = np.where(anyp, DWm[rows, best], np.nan)
            tab["p_localisation"] = np.where(anyp, Pm[rows, best], np.nan)
            tab["localised"] = tab["p_localisation"] <= params.localisation_p
            hit = hit & tab["localised"]
        tab["hit"] = hit.astype(bool)
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
