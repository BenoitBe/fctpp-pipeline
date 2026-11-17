"""Shared condition × temperature component: a diagnostic for effects common to
many proteins.

Why
---
A ligand stabilises one protein or a few. When a large part of the proteome moves
together between conditions at some temperatures, the cause is usually not a
binding event: sample handling per temperature (heating, centrifugation,
pipetting of one condition), a solvent or matrix effect, or the co-aggregation of
large complexes (ribosome). Within-lysate permutations cannot see such an effect,
because it is the same in every lysate: it *is* a condition effect. Only its
shape across temperatures and its lack of specificity tell it apart from a
thermal shift.

Statistic
---------
For a contrast, Δ[p, T] is the protein × temperature matrix of log2 fold changes
(the refits of the shift test). On the proteins quantified at every temperature,
the first right singular vector ``v`` (unit norm, uncentred SVD: a global shift
at one temperature counts) is the shared temperature pattern, and

    variance explained = s1² / Σ sk²

Null: Δ[p, T] → ±Δ[p, T] with independent random signs (each cell keeps its
magnitude, so the heteroscedasticity across temperatures and proteins is kept,
and any structure shared between proteins is destroyed). The component is
*significant* when the observed share exceeds the sign-flip draws
(p = (1 + #null ≥ obs) / (1 + B)), and *strong* when it is also at least
``min_ratio`` (default 1.5) times the 95 % quantile of the draws. With thousands
of proteins a weak structure is always significant; a single noisy temperature
(often the hottest one) can also carry the first component under the null, which
raises the null quantile. Only a strong component flags proteins.

Per protein, ``r_shared`` is the Pearson correlation of its Δ profile with ``v``
(over the temperatures where it is quantified, at least half of the series).
When the component is strong, proteins with |r_shared| ≥ ``r_threshold``
(default 0.7) are flagged ``shared_profile``: their candidate shift follows the
proteome-wide pattern and is not specific. Hits are flagged, never removed.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


@dataclass
class SharedParams:
    enabled: bool = True
    r_threshold: float = 0.7
    n_null: int = 100
    alpha: float = 0.05
    min_ratio: float = 1.5
    seed: int = 42

    def validate(self) -> None:
        if not 0 < self.r_threshold <= 1:
            raise ValueError("shared_component.r_threshold must be in (0, 1]")
        if self.n_null < 10:
            raise ValueError("shared_component.n_null must be >= 10")
        if self.min_ratio < 1:
            raise ValueError("shared_component.min_ratio must be >= 1")


@dataclass
class SharedComponent:
    contrast: str
    temperatures: list[str]
    loadings: pd.Series                 # per temperature, unit norm, largest |loading| positive
    variance_explained: float
    null_q95: float
    p_value: float
    n_proteins: int
    significant: bool
    r: pd.Series                        # per protein (signed)
    second_variance_explained: float = np.nan
    strong: bool = False
    notes: list[str] = field(default_factory=list)

    def flags(self, r_threshold: float) -> pd.Series:
        return (self.r.abs() >= r_threshold) & bool(self.strong)

    @property
    def status(self) -> str:
        if self.strong:
            return "strong"
        return "weak (significant, below min_ratio x null)" if self.significant else "not significant"


def _first_share(M: np.ndarray) -> tuple[float, np.ndarray, float]:
    _, s, vt = np.linalg.svd(M, full_matrices=False)
    tot = float(np.sum(s ** 2))
    if tot <= 0:
        return 0.0, vt[0], 0.0
    return float(s[0] ** 2 / tot), vt[0], float(s[1] ** 2 / tot) if len(s) > 1 else np.nan


def profile_correlation(D: pd.DataFrame, v: pd.Series, min_points: int) -> pd.Series:
    """Pearson r of each row of D with v over the non-missing temperatures."""
    X = D[v.index].to_numpy(float)
    out = np.full(X.shape[0], np.nan)
    vv = v.to_numpy(float)
    for i, row in enumerate(X):
        ok = np.isfinite(row)
        if ok.sum() < min_points:
            continue
        a, b = row[ok] - row[ok].mean(), vv[ok] - vv[ok].mean()
        den = np.sqrt((a ** 2).sum() * (b ** 2).sum())
        if den > 0:
            out[i] = float((a * b).sum() / den)
    return pd.Series(out, index=D.index)


def shared_component(D: pd.DataFrame, contrast: str, params: SharedParams | None = None) -> SharedComponent | None:
    """``D``: proteins × temperatures matrix of log2 FC (NaN where not quantified)."""
    params = params or SharedParams()
    params.validate()
    temps = list(D.columns)
    K = len(temps)
    complete = D.dropna()
    if K < 4 or len(complete) < max(50, 5 * K):
        log.info("Shared component %s: not estimated (%d temperatures, %d complete proteins)",
                 contrast, K, len(complete))
        return None
    M = complete.to_numpy(float)
    ve, v, ve2 = _first_share(M)
    if v[np.argmax(np.abs(v))] < 0:
        v = -v
    rng = np.random.default_rng(params.seed)
    null = np.array([_first_share(M * rng.choice([-1.0, 1.0], size=M.shape))[0]
                     for _ in range(params.n_null)])
    p = float((1 + np.sum(null >= ve)) / (1 + len(null)))
    loadings = pd.Series(v, index=temps)
    r = profile_correlation(D, loadings, max(4, int(np.ceil(K / 2))))
    q95 = float(np.quantile(null, 0.95))
    sig = p <= params.alpha
    sc = SharedComponent(contrast, temps, loadings, ve, q95, p, len(complete), sig, r, ve2,
                         strong=bool(sig and ve >= params.min_ratio * q95))
    log.info("Shared component %s: %.1f %% of the Δ variance (sign-flip null 95 %%: %.1f %%, p = %.3g) — %s%s",
             contrast, 100 * ve, 100 * q95, p, sc.status,
             f"; {int(sc.flags(params.r_threshold).sum())} proteins with |r| ≥ {params.r_threshold:g}"
             if sc.strong else "")
    return sc


def delta_matrix(fits: dict, temperatures: list[str], contrast: str, names: list[str]) -> pd.DataFrame:
    """Proteins × temperatures log2 FC from the per-temperature refits."""
    cols = {}
    for t in temperatures:
        f = fits[t]
        j = f.contrasts.index(contrast)
        cols[t] = pd.Series(f.diff[:, j], index=f.proteins).reindex(names)
    return pd.DataFrame(cols, index=pd.Index(names, name="name"))


def component_table(components: dict[str, SharedComponent | None], params: SharedParams,
                    hits: dict[str, set] | None = None) -> pd.DataFrame:
    """One row per contrast: share, null, p, counts, then the loading of every temperature."""
    rows = []
    for c, sc in components.items():
        if sc is None:
            rows.append({"contrast": c, "status": "not estimated"})
            continue
        flag = sc.flags(params.r_threshold)
        row = {"contrast": c, "status": sc.status,
               "variance_explained": round(sc.variance_explained, 4),
               "ratio_to_null_95pct": round(sc.variance_explained / sc.null_q95, 2) if sc.null_q95 > 0 else None,
               "second_component": round(sc.second_variance_explained, 4),
               "null_95pct": round(sc.null_q95, 4), "p_value": round(sc.p_value, 4),
               "proteins_complete": sc.n_proteins,
               f"proteins_abs_r_ge_{params.r_threshold:g}": int((sc.r.abs() >= params.r_threshold).sum()),
               "share_of_proteome": round(float((sc.r.abs() >= params.r_threshold).mean()), 4)}
        if hits is not None:
            h = hits.get(c, set())
            row["hits"] = len(h)
            row["hits_shared_profile"] = int(flag.reindex(list(h)).fillna(False).sum())
        for t, x in sc.loadings.items():
            row[f"loading_{t}"] = round(float(x), 3)
        rows.append(row)
    return pd.DataFrame(rows)
