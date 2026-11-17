"""Core TPP logic — port of ``tidy_tpp()`` (R/tidy_tpp.R) with explicit rules.

Input: a *long* table, one row per (Temperature, protein), with the per-contrast
statistics produced upstream (Deimos): ``<contrast>_diff`` (log2 FC),
``<contrast>_p.val`` and optionally ``<contrast>_p.adj``.

Steps (same as R):
1. Presence statistics per protein: number of temperatures where the protein
   was quantified (``nb_temp``) and longest run of consecutive temperatures
   (``max_consecutive``).
2. Protein filter: ``nb_temp >= min_temperatures`` and
   ``max_consecutive >= min_streak``.
3. Per contrast: ``<c>_logpval = -log10(p)`` and a per-protein call
   ``<c>_threshold`` = "Over" / "Under".

Call rules
----------
``r_legacy`` reproduces R exactly: "Over" if two *successive rows* of the
protein (sorted by temperature) both pass ``|diff| >= fc_cutoff`` and
``-log10(p) >= logp_cutoff``. Successive rows are not necessarily adjacent
temperatures (a gap at T5 makes T4 and T6 "successive") and the sign of the
fold change is not checked.

``strict`` (default) requires ``min_consecutive`` *adjacent* temperatures of the
series, all passing, with the *same sign* of fold change.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd

from .temperatures import TemperatureSeries

log = logging.getLogger(__name__)

STAT_SUFFIXES = ("diff", "p.adj", "p.val")
CALL_RULES = ("strict", "r_legacy")


# ============================================================================
# Parameters
# ============================================================================

@dataclass
class FilterParams:
    """Presence filter. Integers count temperatures (R defaults: 6 and 6, tuned
    for 12 points); a value < 1 is a fraction of the X temperatures of the
    series (e.g. 0.5 → half of them, rounded up)."""
    min_temperatures: float = 6   # R: min_temperatures
    min_streak: float = 6         # R: min_streak

    @staticmethod
    def _count(value: float, n: int) -> int:
        value = float(value)
        if value < 0:
            raise ValueError("filter thresholds must be >= 0")
        return int(np.ceil(value * n)) if 0 < value < 1 else int(value)

    def resolve(self, n_temperatures: int) -> tuple[int, int]:
        """(min_temperatures, min_streak) in number of temperatures, checked against X."""
        mt = self._count(self.min_temperatures, n_temperatures)
        ms = self._count(self.min_streak, n_temperatures)
        if mt > n_temperatures or ms > n_temperatures:
            raise ValueError(
                f"Presence filter impossible with {n_temperatures} temperatures: "
                f"min_temperatures={mt}, min_streak={ms}. Lower them (e.g. --min-temperatures "
                f"{max(1, n_temperatures // 2)} --min-streak {max(1, n_temperatures // 2)}) or give "
                f"fractions (0.5 = half of the series).")
        return mt, ms


@dataclass
class CallParams:
    fc_cutoff: float = 0.6      # |log2 FC| threshold (R default 0.6)
    logp_cutoff: float = 2.0    # -log10(p) threshold (R default 2, i.e. p < 0.01)
    p_column: str = "p.val"     # "p.val" (R behaviour) or "p.adj"
    min_consecutive: int = 2    # R: 2 (hard-coded)
    rule: str = "strict"        # "strict" | "r_legacy"
    abundance_check: bool = True

    def validate(self) -> None:
        if self.rule not in CALL_RULES:
            raise ValueError(f"call rule must be one of {CALL_RULES}, got {self.rule!r}")
        if self.p_column not in ("p.val", "p.adj"):
            raise ValueError(f"p_column must be 'p.val' or 'p.adj', got {self.p_column!r}")
        if self.min_consecutive < 1:
            raise ValueError("min_consecutive must be >= 1")
        if self.fc_cutoff < 0 or self.logp_cutoff < 0:
            raise ValueError("fc_cutoff and logp_cutoff must be >= 0")


@dataclass
class TPPResult:
    data: pd.DataFrame            # R 'TPP_data' (+ annotations), retained proteins only
    stats: pd.DataFrame           # R 'Stats_proteins' (+ temperatures_present, retained)
    calls: pd.DataFrame           # one row per (contrast, protein): call metrics
    hits: pd.DataFrame            # subset of calls with call == "Over"
    summary: pd.DataFrame         # one row per contrast
    contrasts: list[str]
    temperatures: TemperatureSeries
    filter_params: FilterParams
    call_params: CallParams
    n_lost: int = 0
    notes: list[str] = field(default_factory=list)

    def params_dict(self) -> dict:
        return {"filter": asdict(self.filter_params), "call": asdict(self.call_params)}


# ============================================================================
# Helpers
# ============================================================================

def detect_contrasts(columns) -> list[str]:
    """Contrasts present in a table: ``<c>_p.val`` with a matching ``<c>_diff``.

    R used ``grep("_p\\.adj$")``; requiring the diff/p.val pair is stricter and
    ignores Deimos' secondary ``<c>_p.val_limma`` columns.
    """
    cols = [str(c) for c in columns]
    present = set(cols)
    out: list[str] = []
    for c in cols:
        if c.endswith("_p.val"):
            base = c[: -len("_p.val")]
            if base and f"{base}_diff" in present and base not in out:
                out.append(base)
    return out


def longest_run(positions: np.ndarray) -> int:
    """Longest run of consecutive integers in a sorted array (R: rle(diff == 1))."""
    if len(positions) == 0:
        return 0
    best = cur = 1
    for a, b in zip(positions[:-1], positions[1:]):
        cur = cur + 1 if b - a == 1 else 1
        best = max(best, cur)
    return best


def _truthy(v: object) -> bool:
    """True only for real True values (None/NaN/False -> False)."""
    if v is None:
        return False
    try:
        if isinstance(v, float) and np.isnan(v):
            return False
    except TypeError:
        pass
    return bool(v)


def _neglog10(p: pd.Series) -> pd.Series:
    x = pd.to_numeric(p, errors="coerce").astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return -np.log10(x)


def _segments(mask: np.ndarray, link: np.ndarray) -> list[tuple[int, int]]:
    """Maximal [start, end] index segments where mask is True and consecutive
    elements are linked (link[i] tells whether i and i+1 may be chained)."""
    segs: list[tuple[int, int]] = []
    start = None
    for i, m in enumerate(mask):
        if m and start is None:
            start = i
        if start is not None:
            last = i == len(mask) - 1
            if not m:
                segs.append((start, i - 1))
                start = None
            elif last or not (mask[i + 1] and link[i]):
                segs.append((start, i))
                start = None
    return segs


def _best_segment(segs: list[tuple[int, int]], logp: np.ndarray) -> tuple[int, int] | None:
    if not segs:
        return None
    # longest first, then strongest (max -log10 p), then earliest
    return max(segs, key=lambda s: (s[1] - s[0] + 1, np.nanmax(logp[s[0]:s[1] + 1]), -s[0]))


# ============================================================================
# Steps
# ============================================================================

def presence_stats(data: pd.DataFrame, temps: TemperatureSeries) -> pd.DataFrame:
    """Per-protein ``nb_temp`` and ``max_consecutive`` (R ``stats_prot``)."""
    pos = data["Temperature"].astype(str).map(temps.position)
    df = pd.DataFrame({"name": data["name"].astype(str).to_numpy(),
                       "pos": pos.to_numpy()}).drop_duplicates()
    rows = []
    for name, sub in df.groupby("name", sort=True):
        p = np.sort(sub["pos"].to_numpy(dtype=int))
        rows.append((name, len(p), longest_run(p),
                     ",".join(temps.ids[i] for i in p)))
    return pd.DataFrame(rows, columns=["name", "nb_temp", "max_consecutive",
                                       "temperatures_present"])


def call_protein(pos: np.ndarray, diff: np.ndarray, logp: np.ndarray,
                 params: CallParams, num_nas: np.ndarray | None = None) -> dict:
    """Call one protein for one contrast. Arrays are sorted by temperature."""
    with np.errstate(invalid="ignore"):
        ok = (np.abs(diff) >= params.fc_cutoff) & (logp >= params.logp_cutoff)
    ok = np.where(np.isnan(diff) | np.isnan(logp), False, ok)
    sign = np.sign(np.nan_to_num(diff))
    k = params.min_consecutive
    n = len(pos)

    adjacent = np.diff(pos) == 1 if n > 1 else np.zeros(0, bool)
    same_sign = sign[:-1] == sign[1:] if n > 1 else np.zeros(0, bool)
    link_strict = adjacent & same_sign
    link_legacy = np.ones(max(n - 1, 0), dtype=bool)   # successive rows

    segs_strict = _segments(ok, link_strict)
    segs_legacy = _segments(ok, link_legacy)
    over_strict = any(e - s + 1 >= k for s, e in segs_strict)
    over_legacy = any(e - s + 1 >= k for s, e in segs_legacy)

    segs = segs_strict if params.rule == "strict" else segs_legacy
    over = over_strict if params.rule == "strict" else over_legacy
    best = _best_segment(segs, logp)

    out = {
        "call": "Over" if over else "Under",
        "n_sig_temps": int(ok.sum()),
        "run_length": 0, "run_start": None, "run_end": None,
        "direction": None, "max_abs_diff": np.nan, "max_logp": np.nan,
        "over_strict": bool(over_strict), "over_r_legacy": bool(over_legacy),
        "imputed_in_run": None,
    }
    if best is not None:
        s, e = best
        d = diff[s:e + 1]
        out.update({
            "run_length": int(e - s + 1),
            "run_start": int(pos[s]), "run_end": int(pos[e]),
            "direction": "up" if np.all(d > 0) else ("down" if np.all(d < 0) else "mixed"),
            "max_abs_diff": float(np.nanmax(np.abs(d))),
            "max_logp": float(np.nanmax(logp[s:e + 1])),
        })
        if num_nas is not None:
            nn = num_nas[s:e + 1]
            out["imputed_in_run"] = bool(np.nansum(nn) > 0) if not np.all(np.isnan(nn)) else None
    return out


def run_tpp(results: pd.DataFrame, temps: TemperatureSeries | None = None,
            filter_params: FilterParams | None = None,
            call_params: CallParams | None = None,
            contrasts: list[str] | None = None) -> TPPResult:
    """Full TPP layer on a long results table (see module docstring)."""
    fp = filter_params or FilterParams()
    cp = call_params or CallParams()
    cp.validate()

    if not {"Temperature", "name"}.issubset(results.columns):
        raise ValueError("results must contain 'Temperature' and 'name' columns")
    data = results.copy()
    data["Temperature"] = data["Temperature"].astype(str)
    data["name"] = data["name"].astype(str)
    if temps is None:
        temps = TemperatureSeries.from_ids(data["Temperature"].unique())
    unknown = sorted(set(data["Temperature"]) - set(temps.ids))
    if unknown:
        raise ValueError(f"Temperatures {unknown} are not in the series {temps.ids}")
    dup = data.duplicated(["Temperature", "name"])
    if dup.any():
        ex = data.loc[dup, ["Temperature", "name"]].head(5).to_dict("records")
        raise ValueError(f"{int(dup.sum())} duplicated (Temperature, name) rows, e.g. {ex}")

    contrasts = contrasts or detect_contrasts(data.columns)
    if not contrasts:
        raise ValueError("No contrast found (expected '<contrast>_diff' + '<contrast>_p.val' columns)")
    if cp.p_column == "p.adj":
        missing = [c for c in contrasts if f"{c}_p.adj" not in data.columns]
        if missing:
            raise ValueError(f"p_column='p.adj' but no '_p.adj' column for {missing}")

    data["_pos"] = data["Temperature"].map(temps.position).astype(int)
    data = data.sort_values(["_pos", "name"], kind="mergesort").reset_index(drop=True)

    # --- 1-2. presence statistics & filter -----------------------------------
    min_t, min_s = fp.resolve(len(temps))
    if cp.min_consecutive > len(temps):
        raise ValueError(f"min_consecutive={cp.min_consecutive} > {len(temps)} temperatures")
    stats = presence_stats(data, temps)
    keep = (stats["max_consecutive"] >= min_s) & (stats["nb_temp"] >= min_t)
    stats["retained"] = keep.to_numpy()
    kept_names = set(stats.loc[keep, "name"])
    tpp = data[data["name"].isin(kept_names)].copy()
    n_lost = data["name"].nunique() - tpp["name"].nunique()
    log.info("%d temperatures; %d proteins lost after the presence filter (nb_temp >= %d, "
             "streak >= %d); %d kept.", len(temps), n_lost, min_t, min_s, len(kept_names))
    if not kept_names:
        log.warning("No protein passes the presence filter — check filter.min_temperatures / min_streak")

    # --- 3. calls -------------------------------------------------------------
    has_nas = "num_NAs" in tpp.columns
    calls_rows = []
    for c in contrasts:
        tpp[f"{c}_logpval"] = _neglog10(tpp[f"{c}_{cp.p_column}"])
        tpp[f"{c}_threshold"] = "Under"

    for name, sub in tpp.groupby("name", sort=True):
        sub = sub.sort_values("_pos")
        pos = sub["_pos"].to_numpy(dtype=int)
        nn = (pd.to_numeric(sub["num_NAs"], errors="coerce").to_numpy(dtype=float)
              if has_nas else None)
        for c in contrasts:
            diff = pd.to_numeric(sub[f"{c}_diff"], errors="coerce").to_numpy(dtype=float)
            logp = sub[f"{c}_logpval"].to_numpy(dtype=float)
            res = call_protein(pos, diff, logp, cp, nn)
            res.update({"contrast": c, "name": name})
            # abundance check at the reference (lowest) temperature
            res["reference_temperature"] = temps.reference
            res["shift_class"] = None
            if res["call"] == "Over" and cp.abundance_check:
                ref = np.where(pos == 0)[0]
                if res["direction"] == "mixed":
                    res["shift_class"] = "mixed-sign"
                elif len(ref) == 0:
                    res["shift_class"] = "undetermined (absent at reference)"
                else:
                    i = ref[0]
                    ref_ok = (abs(diff[i]) >= cp.fc_cutoff) and (logp[i] >= cp.logp_cutoff)
                    same = res["direction"] in ("up", "down") and \
                        np.sign(diff[i]) == (1 if res["direction"] == "up" else -1)
                    res["shift_class"] = "abundance-like" if (ref_ok and same) else "thermal"
            calls_rows.append(res)

    calls = pd.DataFrame(calls_rows)
    if len(calls):
        for col in ("run_start", "run_end"):
            calls[col] = calls[col].map(lambda p: temps.ids[int(p)] if pd.notna(p) else None)
        over = calls.loc[calls["call"] == "Over", ["contrast", "name"]]
        for c in contrasts:
            names_over = set(over.loc[over["contrast"] == c, "name"])
            tpp.loc[tpp["name"].isin(names_over), f"{c}_threshold"] = "Over"

    # annotations for readability
    ann_cols = [a for a in ("Genes", "Protein.Group", "First.Protein.Description") if a in tpp.columns]
    if ann_cols and len(calls):
        ann = (tpp.sort_values("_pos").groupby("name")[ann_cols].first().reset_index())
        calls = calls.merge(ann, on="name", how="left")

    front = ["contrast", "name"] + ann_cols + [
        "call", "shift_class", "direction", "run_start", "run_end", "run_length",
        "n_sig_temps", "max_abs_diff", "max_logp", "over_strict", "over_r_legacy",
        "imputed_in_run", "reference_temperature"]
    calls = calls[[c for c in front if c in calls.columns]] if len(calls) else pd.DataFrame(columns=front)
    hits = calls[calls["call"] == "Over"].reset_index(drop=True)

    # --- R-compatible TPP_data layout ----------------------------------------
    stat_cols = []
    for c in contrasts:
        for s in STAT_SUFFIXES:
            if f"{c}_{s}" in tpp.columns:
                stat_cols.append(f"{c}_{s}")
    derived = [f"{c}_{s}" for c in contrasts for s in ("logpval", "threshold")]
    extra = [x for x in ("Temperature_C",) if x in tpp.columns]
    ann_tail = [x for x in ("Genes", "First.Protein.Description", "num_NAs", "imputed")
                if x in tpp.columns]
    tpp = tpp[["Temperature"] + extra + ["name"] + stat_cols + derived + ann_tail].reset_index(drop=True)

    summary = _summarise(calls, contrasts, len(kept_names), has_nas)
    for c in contrasts:
        stars = hits.loc[hits["contrast"] == c, "name"].tolist()
        log.info("Contrast '%s': %d proteins passed the cutoff%s", c, len(stars),
                 (": " + ", ".join(stars[:30]) + (" ..." if len(stars) > 30 else "")) if stars else ".")

    return TPPResult(data=tpp, stats=stats, calls=calls, hits=hits, summary=summary,
                     contrasts=list(contrasts), temperatures=temps,
                     filter_params=fp, call_params=cp, n_lost=int(n_lost))


def _summarise(calls: pd.DataFrame, contrasts: list[str], n_kept: int,
               has_nas: bool) -> pd.DataFrame:
    rows = []
    for c in contrasts:
        sub = calls[calls["contrast"] == c] if len(calls) else calls
        h = sub[sub["call"] == "Over"] if len(sub) else sub
        rows.append({
            "contrast": c,
            "proteins_tested": int(len(sub)),
            "hits": int(len(h)),
            "hits_up": int((h["direction"] == "up").sum()) if len(h) else 0,
            "hits_down": int((h["direction"] == "down").sum()) if len(h) else 0,
            "hits_thermal": int((h["shift_class"] == "thermal").sum()) if len(h) else 0,
            "hits_abundance_like": int((h["shift_class"] == "abundance-like").sum()) if len(h) else 0,
            "hits_with_imputed_values": ((int(sum(_truthy(v) for v in h["imputed_in_run"])) if len(h) else 0)
                                         if has_nas else "n/a"),
            "hits_strict_rule": int(sub["over_strict"].sum()) if len(sub) else 0,
            "hits_r_legacy_rule": int(sub["over_r_legacy"].sum()) if len(sub) else 0,
        })
    out = pd.DataFrame(rows)
    out.attrs["proteins_retained"] = n_kept
    return out


def wide_matrix(tpp: TPPResult, contrast: str, value: str = "diff") -> pd.DataFrame:
    """Protein × temperature matrix of ``<contrast>_<value>`` (diff, logpval, p.val…)."""
    col = f"{contrast}_{value}"
    m = tpp.data.pivot(index="name", columns="Temperature", values=col)
    return m.reindex(columns=[t for t in tpp.temperatures.ids if t in m.columns])
