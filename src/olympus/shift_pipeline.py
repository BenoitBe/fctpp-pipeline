"""Shift test orchestration: matrices → paired refits → protein-level test →
permutation FDR estimates. Called by ``pipeline.analyze``."""

from __future__ import annotations

import logging
from contextlib import contextmanager
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .engine import fit_series, load_limma
from .matrices import has_matrices, load_series_matrices
from .shift import (ShiftParams, adjacent_correlation, fits_to_long, permutation_null,
                    shift_test, swap_patterns)
from .tpp import CallParams, TPPResult, run_tpp

log = logging.getLogger(__name__)


@dataclass
class ShiftResult:
    table: pd.DataFrame                     # protein level, all contrasts
    windows: pd.DataFrame                   # window level
    fdr: pd.DataFrame                       # permutation estimates per contrast × rule
    params: ShiftParams
    paired_used: bool
    rho: dict
    fits: dict
    legacy_paired: bool
    n_patterns: int
    notes: list[str] = field(default_factory=list)


@contextmanager
def _quiet(*names):
    loggers = [logging.getLogger(n) for n in names]
    old = [lg.level for lg in loggers]
    for lg in loggers:
        lg.setLevel(logging.WARNING)
    try:
        yield
    finally:
        for lg, lv in zip(loggers, old):
            lg.setLevel(lv)


def _deimos_design_paired(methods: pd.DataFrame | None) -> bool:
    """True when Deimos fitted ~0 + condition + subject (Methods 'Design matrix')."""
    if methods is None or not len(methods) or "item" not in methods.columns:
        return False
    row = methods[methods["item"] == "Design matrix"]
    if row.empty:
        return False
    vals = [str(v) for v in row.iloc[0].drop("item").dropna()]
    return any("subject" in v for v in vals)


def _legacy_with_pvalue(cp: CallParams) -> CallParams:
    return CallParams(fc_cutoff=cp.fc_cutoff, logp_cutoff=cp.logp_cutoff, p_column=cp.p_column,
                      min_consecutive=cp.min_consecutive, rule=cp.rule, abundance_check=False)


def run_shift_analysis(cfg: dict, loaded, temps, res: TPPResult, mode: str) -> ShiftResult | None:
    st = cfg.get("shift_test", {}) or {}
    enabled = st.get("enabled", "auto")
    if enabled in (False, "false", "False", "off"):
        return None
    strict = enabled in (True, "true", "True", "on")

    def _skip(msg):
        if strict:
            raise RuntimeError(msg)
        log.info("Shift test skipped: %s", msg)
        return None

    if mode != "deimos":
        return _skip("needs the Deimos workbooks (input.mode: deimos)")
    runs = dict(zip(loaded.provenance["Temperature"].astype(str), loaded.provenance["path"]))
    missing = [t for t in temps.ids if not has_matrices(runs.get(t, ""))]
    if missing:
        return _skip(f"no Deimos matrices/ folder for {missing} (export_matrices: true in Deimos)")
    try:
        limma = load_limma(st.get("deimos_dir"))
    except ImportError as e:
        return _skip(str(e))

    cp = res.call_params
    params = ShiftParams(
        baseline=str(st.get("baseline", "none")),
        window=int(st.get("window") or cp.min_consecutive),
        alpha=float(st.get("alpha", 0.05)),
        fc_cutoff=float(st.get("fc_cutoff") if st.get("fc_cutoff") is not None else cp.fc_cutoff),
        use_adjacent_correlation=bool(st.get("use_adjacent_correlation", True)),
        paired=st.get("paired", "auto"),
        localisation_p=(None if st.get("localisation_p", 0.01) in (None, False, "null")
                        else float(st.get("localisation_p", 0.01))),
    )
    params.validate()
    mats = load_series_matrices({t: runs[t] for t in temps.ids})
    universe = sorted(res.stats.loc[res.stats["retained"], "name"].astype(str))
    uni = set(universe)
    contrasts = list(res.contrasts)

    fits = fit_series(limma, mats, list(temps.ids), contrasts, params.paired, proteins=uni)
    paired_used = all(f.paired for f in fits.values())
    rho = adjacent_correlation(fits, list(temps.ids))
    table, windows = shift_test(fits, list(temps.ids), contrasts, params, rho, universe)
    ann = res.data.groupby("name")[[c for c in ("Genes", "First.Protein.Description")
                                    if c in res.data.columns]].first()
    if len(ann.columns):
        table = table.merge(ann, left_on="name", right_index=True, how="left")
    log.info("Shift test (%s fit, baseline %s, window %d): %s", "paired" if paired_used else "unpaired",
             params.baseline, params.window,
             ", ".join(f"{c}: {int(table[(table.contrast == c)].hit.sum())} hits" for c in contrasts))

    legacy_paired = _deimos_design_paired(getattr(loaded, "methods", None))
    fdr_rows, n_patterns = [], 0
    if st.get("permutations", True):
        blocks = sorted(set(mats[temps.ids[0]].block_labels()))
        n_patterns = len(swap_patterns(blocks, int(st.get("max_patterns", 64)), int(st.get("seed", 42))))
        leg_cp = _legacy_with_pvalue(cp)

        def eval_shift(f, c):
            t2, _ = shift_test(f, list(temps.ids), [c], params, rho, universe)
            return int(t2["hit"].sum())

        def eval_legacy(f, c):
            long = fits_to_long(f, list(temps.ids), [c])
            r = run_tpp(long, temps, res.filter_params, leg_cp, contrasts=[c])
            return int(len(r.hits))

        kw = dict(max_patterns=int(st.get("max_patterns", 64)), seed=int(st.get("seed", 42)), proteins=uni)
        with _quiet("olympus.tpp"):
            null_shift = permutation_null(limma, mats, list(temps.ids), contrasts, eval_shift,
                                          params.paired, **kw)
            fits_leg = fit_series(limma, mats, list(temps.ids), contrasts, legacy_paired, proteins=uni)
            obs_leg = {c: eval_legacy(fits_leg, c) for c in contrasts}
            null_leg = permutation_null(limma, mats, list(temps.ids), contrasts, eval_legacy,
                                        legacy_paired, **kw)
        for c in contrasts:
            for rule, obs, null, fit in (
                    ("shift test (q ≤ %.2g)" % params.alpha, int(table[table.contrast == c].hit.sum()),
                     null_shift[c], "paired" if paired_used else "unpaired"),
                    (f"legacy rule ({cp.rule})", obs_leg[c], null_leg[c], "paired" if legacy_paired else "unpaired")):
                m = float(np.mean(null)) if null else np.nan
                fdr_rows.append({
                    "contrast": c, "rule": rule, "fit": fit, "observed_hits": obs,
                    "null_hits_mean": round(m, 2) if null else None,
                    "null_hits": ", ".join(map(str, null)),
                    "fdr_estimate": (round(min(1.0, m / obs), 2) if obs > 0 else (None if not null else 0.0)),
                    "patterns": len(null)})
        legacy_from_deimos = {c: int((res.hits["contrast"] == c).sum()) for c in contrasts}
        for c in contrasts:
            if legacy_from_deimos[c] != obs_leg[c]:
                log.info("Legacy rule %s: %d hits from the Deimos sheet, %d from the refit "
                         "(rounding of matrices/ or different fit)", c, legacy_from_deimos[c], obs_leg[c])
    fdr = pd.DataFrame(fdr_rows)
    return ShiftResult(table.sort_values(["contrast", "p_protein"], kind="mergesort").reset_index(drop=True),
                       windows, fdr, params, paired_used, rho, fits, legacy_paired, n_patterns)


def summary_columns(sr: ShiftResult | None, res: TPPResult) -> pd.DataFrame:
    """Columns added to the Summary sheet."""
    s = res.summary.copy()
    if sr is None:
        return s
    s["shift_hits"] = s["contrast"].map(sr.table.groupby("contrast")["hit"].sum()).fillna(0).astype(int)
    s["shift_min_q"] = s["contrast"].map(sr.table.groupby("contrast")["q_protein"].min()).round(4)
    if len(sr.fdr):
        leg = sr.fdr[sr.fdr["rule"].str.startswith("legacy")].set_index("contrast")
        s["legacy_rule_fdr_estimate"] = s["contrast"].map(leg["fdr_estimate"])
    return s


def methods_rows(sr: ShiftResult) -> list[tuple[str, str, str]]:
    p = sr.params
    rho_txt = ", ".join(f"{a}-{b}: {v:.2f}" for (a, b), v in sr.rho.items())
    return [
        ("Shift test", f"{'paired' if sr.paired_used else 'unpaired'} refit, baseline {p.baseline}, "
                       f"window {p.window}",
         "Deimos limma (limma_ebayes.py) refitted at each temperature on matrices/log2_imputed.tsv"
         + (" with the lysate (replicate) as block: ~0 + condition + replicate" if sr.paired_used else
            ": ~0 + condition")
         + f". Window statistic D = mean(Δ over {p.window} adjacent temperatures) − baseline "
           f"({p.baseline}); variance includes the adjacent-temperature correlation; Welch–Satterthwaite df."),
        ("Shift test FDR", f"Simes over windows, BH over proteins, q ≤ {p.alpha}, |D| ≥ {p.fc_cutoff}",
         "Protein-level false discovery rate, per contrast, over the proteins retained by the presence filter"),
        ("Localisation check", ("off" if p.localisation_p is None or p.baseline == "median"
                                else f"best window vs profile median, p ≤ {p.localisation_p}"),
         "Rejects constant offsets (significant at every temperature, never localised): handling "
         "artefacts in lysate experiments, abundance changes in cells"),
        ("Adjacent correlation", "consensus residual correlation, same lysate", rho_txt),
        ("Permutations", f"{sr.n_patterns} within-lysate swap pattern(s) per contrast",
         "A/B labels swapped in a subset of lysates, the same at every temperature (identity and full "
         "swap excluded; a pattern and its complement are equivalent). FDR estimate = mean null hits / "
         "observed hits. Indicative with few lysates; partial swaps keep part of a true signal (conservative)."),
    ]
