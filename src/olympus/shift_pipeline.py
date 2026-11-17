"""Shift test orchestration: matrices → paired refits → protein-level test →
permutation FDR estimates. Called by ``pipeline.analyze``."""

from __future__ import annotations

import logging
from contextlib import contextmanager
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .engine import BlockCorrelation, block_correlation, fit_series, load_limma
from .matrices import has_matrices, load_series_matrices
from .shared import SharedParams, component_table, delta_matrix, shared_component
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
    block: BlockCorrelation | None = None
    paired_reason: str = ""
    shared: dict = field(default_factory=dict)          # contrast -> SharedComponent | None
    shared_params: SharedParams | None = None
    shared_table: pd.DataFrame | None = None
    notes: list[str] = field(default_factory=list)

    def shared_flags(self, contrast: str) -> pd.Series:
        sc = self.shared.get(contrast)
        if sc is None or self.shared_params is None:
            return pd.Series(dtype=bool)
        return sc.flags(self.shared_params.r_threshold)


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
    win, wins = st.get("window"), st.get("windows")
    if isinstance(wins, (int, float)):
        wins = [int(wins)]
    if win:                                   # explicit single size wins
        win_kw = {"window": int(win)}
    elif wins:
        win_kw = {"window": None, "windows": tuple(int(w) for w in wins)}
    else:                                     # legacy behaviour: call.min_consecutive
        win_kw = {"window": int(cp.min_consecutive)}
    params = ShiftParams(
        baseline=str(st.get("baseline", "none")),
        **win_kw,
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

    # pairing: forced, or decided from the measured within-lysate correlation
    block = block_correlation(mats, list(temps.ids), uni)
    pv = params.paired
    if pv in (True, "true", "True"):
        paired_flag, reason = True, "forced (paired: true); " + block.describe()
        if not block.complete:
            raise ValueError("shift_test.paired: true but the lysate blocks are incomplete")
    elif pv in (False, "false", "False"):
        paired_flag, reason = False, "forced (paired: false); " + block.describe()
    else:
        paired_flag = block.decide()
        reason = ("auto: " + ("paired" if paired_flag else "unpaired") + " — " + block.describe())
    log.info("Pairing %s", reason)

    fits = fit_series(limma, mats, list(temps.ids), contrasts, paired_flag, proteins=uni)
    paired_used = all(f.paired for f in fits.values())
    rho = adjacent_correlation(fits, list(temps.ids))
    table, windows = shift_test(fits, list(temps.ids), contrasts, params, rho, universe)

    # shared condition x temperature component (diagnostic, flags hits)
    shp = SharedParams(enabled=bool(st.get("shared_component", True)),
                       r_threshold=float(st.get("shared_r_threshold", 0.7)),
                       n_null=int(st.get("shared_n_null", 100)), min_ratio=float(st.get("shared_min_ratio", 1.5)),
                       seed=int(st.get("seed", 42)))
    shared, shared_tab = {}, None
    if shp.enabled:
        shp.validate()
        for c in contrasts:
            shared[c] = shared_component(delta_matrix(fits, list(temps.ids), c, universe), c, shp)
        r_all = pd.concat([sc.r.rename("r_shared").to_frame().assign(contrast=c)
                           for c, sc in shared.items() if sc is not None]) if any(shared.values()) else None
        if r_all is not None:
            table = table.merge(r_all.reset_index(), on=["contrast", "name"], how="left")
            sig = table["contrast"].map({c: bool(sc is not None and sc.strong) for c, sc in shared.items()})
            table["shared_profile"] = (table["r_shared"].abs() >= shp.r_threshold) & sig.fillna(False).astype(bool)
        else:
            table["r_shared"], table["shared_profile"] = np.nan, False
        table["hit_specific"] = table["hit"] & ~table["shared_profile"]
        shared_tab = component_table(shared, shp, {c: set(table.loc[(table.contrast == c) & table.hit, "name"])
                                                   for c in contrasts})
    ann = res.data.groupby("name")[[c for c in ("Genes", "First.Protein.Description")
                                    if c in res.data.columns]].first()
    if len(ann.columns):
        table = table.merge(ann, left_on="name", right_index=True, how="left")
    log.info("Shift test (%s fit, baseline %s, %s): %s", "paired" if paired_used else "unpaired",
             params.baseline, params.describe_windows(),
             ", ".join(f"{c}: {int(table[(table.contrast == c)].hit.sum())} hits"
                       + (f" ({int(table[(table.contrast == c)].hit_specific.sum())} not following the "
                          f"shared component)" if "hit_specific" in table.columns else "")
                       for c in contrasts))

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
                                          paired_flag, **kw)
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
                       windows, fdr, params, paired_used, rho, fits, legacy_paired, n_patterns,
                       block=block, paired_reason=reason, shared=shared, shared_params=shp if shp.enabled else None,
                       shared_table=shared_tab)


def annotate_legacy_hits(sr: ShiftResult | None, hits: pd.DataFrame) -> pd.DataFrame:
    """Add r_shared / shared_profile to the legacy-rule hits (same refits as the shift test)."""
    if sr is None or not sr.shared or hits is None or hits.empty:
        return hits
    out = hits.copy()
    r, flag = [], []
    for c, n in zip(out["contrast"], out["name"]):
        sc = sr.shared.get(c)
        rv = float(sc.r.get(n, np.nan)) if sc is not None else np.nan
        r.append(rv)
        flag.append(bool(sc is not None and sc.strong and np.isfinite(rv)
                         and abs(rv) >= sr.shared_params.r_threshold))
    out["r_shared"] = np.round(r, 3)
    out["shared_profile"] = flag
    return out


def summary_columns(sr: ShiftResult | None, res: TPPResult) -> pd.DataFrame:
    """Columns added to the Summary sheet."""
    s = res.summary.copy()
    if sr is None:
        return s
    s["shift_hits"] = s["contrast"].map(sr.table.groupby("contrast")["hit"].sum()).fillna(0).astype(int)
    if "hit_specific" in sr.table.columns:
        s["shift_hits_not_shared"] = s["contrast"].map(
            sr.table.groupby("contrast")["hit_specific"].sum()).fillna(0).astype(int)
    s["shift_min_q"] = s["contrast"].map(sr.table.groupby("contrast")["q_protein"].min()).round(4)
    if sr.shared:
        s["shared_component_share"] = s["contrast"].map(
            {c: (round(sc.variance_explained, 3) if sc is not None else None) for c, sc in sr.shared.items()})
        s["shared_component_p"] = s["contrast"].map(
            {c: (round(sc.p_value, 3) if sc is not None else None) for c, sc in sr.shared.items()})
        lh = annotate_legacy_hits(sr, res.hits)
        if lh is not None and "shared_profile" in lh.columns:
            s["hits_shared_profile"] = s["contrast"].map(lh.groupby("contrast")["shared_profile"].sum()) \
                .fillna(0).astype(int)
    if len(sr.fdr):
        leg = sr.fdr[sr.fdr["rule"].str.startswith("legacy")].set_index("contrast")
        s["legacy_rule_fdr_estimate"] = s["contrast"].map(leg["fdr_estimate"])
    return s


def methods_rows(sr: ShiftResult) -> list[tuple[str, str, str]]:
    p = sr.params
    rho_txt = ", ".join(f"{a}-{b}: {v:.2f}" for (a, b), v in sr.rho.items())
    rows = [
        ("Pairing", "paired (lysate as block)" if sr.paired_used else "unpaired", sr.paired_reason),
        ("Shift test", f"{'paired' if sr.paired_used else 'unpaired'} refit, baseline {p.baseline}, "
                       f"{p.describe_windows()}",
         "Deimos limma (limma_ebayes.py) refitted at each temperature on matrices/log2_imputed.tsv"
         + (" with the lysate (replicate) as block: ~0 + condition + replicate" if sr.paired_used else
            ": ~0 + condition")
         + f". Window statistic D = mean(Δ over w adjacent temperatures, w ∈ {list(p.window_sizes)}) − "
           f"baseline ({p.baseline}); variance includes the adjacent-temperature correlation; "
           f"Welch–Satterthwaite df."),
        ("Shift test FDR", f"Simes over windows, BH over proteins, q ≤ {p.alpha}, |dd_peak| ≥ {p.fc_cutoff}",
         "Protein-level false discovery rate, per contrast, over the proteins retained by the presence filter"),
        ("Localisation check", ("off" if p.localisation_p is None or p.baseline == "median"
                                else f"best window vs profile median, p ≤ {p.localisation_p}"),
         "Rejects constant offsets (significant at every temperature, never localised): handling "
         "artefacts in lysate experiments, abundance changes in cells"),
        ("Adjacent correlation", "consensus residual correlation, same lysate", rho_txt),
        ("Permutations", f"{sr.n_patterns} within-lysate swap pattern(s) per contrast",
         "A/B labels swapped in a subset of lysates, the same at every temperature (identity and full "
         "swap excluded; a pattern and its complement are equivalent). FDR estimate = mean null hits / "
         "observed hits. Indicative with few lysates; partial swaps keep part of a true signal (conservative). "
         "Blind to effects shared by every lysate: see the shared component."),
    ]
    if sr.shared_params is not None:
        txt = "; ".join(f"{c}: " + ("not estimated" if sc is None else
                                    f"{100 * sc.variance_explained:.1f} % (null 95 % {100 * sc.null_q95:.1f} %, "
                                    f"p = {sc.p_value:.2g}, {sc.status})") for c, sc in sr.shared.items())
        rows.append(("Shared component",
                     f"1st SVD component of the protein × temperature Δ matrix; flag |r| ≥ "
                     f"{sr.shared_params.r_threshold:g} when strong (p ≤ 0.05 and ≥ "
                     f"{sr.shared_params.min_ratio:g}× the null 95 % quantile)",
                     "Effect common to many proteins (handling per temperature, solvent, complexes); "
                     "sign-flip null; hits following it are flagged shared_profile, not removed. " + txt))
    return rows
