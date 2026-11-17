"""Static figures (matplotlib): volcano plots and temperature profiles.

* ``volcano_pdf``  — port of ``magma()``. Default layout: one page per contrast
  with one panel per temperature (shared axes), so a shift can be followed along
  the series; ``layout="pages"`` gives one volcano per page as in R.
  The thresholds default to the TPP call thresholds (in R the volcano used
  FC 0.3 while the calls used FC 0.6).
* ``profiles_pdf`` — port of the 2D part of ``plot_and_3D()``. All retained
  proteins are drawn as a grey background and the hits are highlighted (instead
  of ~1,500 lines each with its own colour and no legend); hits then get one
  small panel each.
"""

from __future__ import annotations

import logging
import math
from contextlib import contextmanager
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from .tpp import TPPResult  # noqa: E402

log = logging.getLogger(__name__)

# Diverging pair (validated: CVD ΔE 24.7, normal ΔE 33.6) + neutral grey.
UP = "#eb6834"
DOWN = "#2a78d6"
NEUTRAL = "#c3c2b7"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"

try:  # optional, ggrepel equivalent
    from adjustText import adjust_text  # type: ignore
except Exception:  # pragma: no cover
    adjust_text = None


@contextmanager
def _style():
    with plt.rc_context({
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
        "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "axes.labelcolor": INK_2,
        "xtick.color": INK_2, "ytick.color": INK_2, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.titleweight": "bold", "text.color": INK, "legend.frameon": False,
        "pdf.fonttype": 42,
    }):
        yield


def _label_series(df: pd.DataFrame, label_column: str | None) -> pd.Series:
    """Gene name when available (more readable than the accession used in R)."""
    lab = df["name"].astype(str)
    if label_column and label_column in df.columns:
        g = df[label_column].astype("string").str.split(";").str[0]
        lab = g.where(g.notna() & (g.str.strip() != ""), lab).astype(str)
    return lab


def _place_labels(ax, xs, ys, texts, max_labels: int = 15):
    if len(texts) == 0:
        return
    order = np.argsort(-np.asarray(ys))[:max_labels]
    objs = [ax.text(xs[i], ys[i], texts[i], fontsize=6, color=INK, ha="center", va="bottom")
            for i in order]
    if adjust_text is not None and len(objs) > 1:
        try:
            adjust_text(objs, ax=ax, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.4))
        except Exception:  # pragma: no cover - cosmetic only
            pass


def _classify(diff: np.ndarray, logp: np.ndarray, fc: float, lp: float) -> np.ndarray:
    cls = np.full(len(diff), "ns", dtype=object)
    with np.errstate(invalid="ignore"):
        cls[(diff >= fc) & (logp >= lp)] = "up"
        cls[(diff <= -fc) & (logp >= lp)] = "down"
    return cls


def _volcano_panel(ax, sub: pd.DataFrame, c: str, fc: float, lp: float, hits: set,
                   label_column: str | None, label_mode: str, max_labels: int):
    diff = pd.to_numeric(sub[f"{c}_diff"], errors="coerce").to_numpy(float)
    logp = pd.to_numeric(sub[f"{c}_logpval"], errors="coerce").to_numpy(float)
    cls = _classify(diff, logp, fc, lp)
    for k, col, z in (("ns", NEUTRAL, 1), ("down", DOWN, 2), ("up", UP, 2)):
        m = cls == k
        ax.scatter(diff[m], logp[m], s=6 if k == "ns" else 9, c=col, linewidths=0,
                   alpha=0.6 if k == "ns" else 0.9, zorder=z, rasterized=True)
    is_hit = sub["name"].isin(hits).to_numpy()
    if is_hit.any():
        ax.scatter(diff[is_hit], logp[is_hit], s=22, facecolors="none", edgecolors=INK,
                   linewidths=0.6, zorder=3)
    for v in (-fc, fc):
        ax.axvline(v, color=MUTED, lw=0.7, ls="--", zorder=1.6)
    ax.axhline(lp, color=MUTED, lw=0.7, ls="--", zorder=1.6)
    labels = _label_series(sub, label_column).to_numpy()
    sel = (is_hit & (cls != "ns")) if label_mode == "hits" else (cls != "ns")
    idx = np.where(sel)[0]
    _place_labels(ax, diff[idx], logp[idx], labels[idx], max_labels)


def volcano_pdf(res: TPPResult, out_pdf: str | Path, *, fc_cutoff: float | None = None,
                logp_cutoff: float | None = None, layout: str = "grid",
                label_column: str | None = "Genes", source: pd.DataFrame | None = None) -> Path:
    """Volcano plots per temperature and contrast.

    ``source`` defaults to ``res.data`` (retained proteins); pass the full long
    table to draw every quantified protein, as R did.
    """
    fc = res.call_params.fc_cutoff if fc_cutoff is None else fc_cutoff
    lp = res.call_params.logp_cutoff if logp_cutoff is None else logp_cutoff
    data = res.data if source is None else source
    temps = res.temperatures
    out_pdf = Path(out_pdf)
    legend = [Line2D([], [], ls="", marker="o", ms=4, color=UP, label="Up"),
              Line2D([], [], ls="", marker="o", ms=4, color=DOWN, label="Down"),
              Line2D([], [], ls="", marker="o", ms=4, color=NEUTRAL, label="Not significant"),
              Line2D([], [], ls="", marker="o", ms=5, mfc="none", mec=INK, label="TPP hit")]
    with _style(), PdfPages(out_pdf) as pdf:
        for c in res.contrasts:
            if f"{c}_logpval" not in data.columns:
                data = data.assign(**{f"{c}_logpval": -np.log10(
                    pd.to_numeric(data[f"{c}_{res.call_params.p_column}"], errors="coerce"))})
            hits = set(res.hits.loc[res.hits["contrast"] == c, "name"])
            present = [t for t in temps.ids if (data["Temperature"] == t).any()]
            if layout == "pages":
                for t in present:
                    fig, ax = plt.subplots(figsize=(7, 5))
                    sub = data[data["Temperature"] == t]
                    _volcano_panel(ax, sub, c, fc, lp, hits, label_column, "all", 15)
                    ax.set_title(f"{temps.label(t)} — {c}")
                    ax.set_xlabel("log2 fold change")
                    ax.set_ylabel(f"-log10({res.call_params.p_column})")
                    fig.legend(handles=legend, loc="upper right", fontsize=7)
                    pdf.savefig(fig)
                    plt.close(fig)
                continue
            n = len(present)
            ncol = min(4, n)
            nrow = math.ceil(n / ncol)
            fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol + 0.6, 2.3 * nrow + 0.9),
                                     sharex=True, sharey=True, squeeze=False)
            for ax, t in zip(axes.flat, present):
                sub = data[data["Temperature"] == t]
                _volcano_panel(ax, sub, c, fc, lp, hits, label_column, "hits", 8)
                ax.set_title(temps.label(t), fontsize=8)
            for ax in list(axes.flat)[n:]:
                ax.set_visible(False)
            for ax in axes[-1, :]:
                ax.set_xlabel("log2 fold change")
            for ax in axes[:, 0]:
                ax.set_ylabel(f"-log10({res.call_params.p_column})")
            fig.suptitle(f"{c} — |log2FC| ≥ {fc:g}, -log10(p) ≥ {lp:g}  ({len(hits)} TPP hits)",
                         fontsize=10, fontweight="bold", color=INK)
            fig.legend(handles=legend, loc="upper center", ncol=4, fontsize=7,
                       bbox_to_anchor=(0.5, 0.965))
            fig.tight_layout(rect=(0, 0, 1, 0.925))
            pdf.savefig(fig, dpi=200)
            plt.close(fig)
    log.info("Volcano plots: %s", out_pdf)
    return out_pdf


def _series_arrays(sub: pd.DataFrame, col: str, temps) -> tuple[np.ndarray, np.ndarray]:
    sub = sub.sort_values("Temperature", key=lambda s: s.map(temps.position))
    x = np.array([temps.x(t) for t in sub["Temperature"]])
    y = pd.to_numeric(sub[col], errors="coerce").to_numpy(float)
    return x, y


def _set_temperature_axis(ax, temps):
    xs = [temps.x(t) for t in temps.ids]
    ax.set_xticks(xs)
    ax.set_xticklabels(temps.tick_labels(), rotation=0 if temps.celsius else 45, fontsize=7)
    ax.set_xlabel(temps.axis_title)


def profiles_pdf(res: TPPResult, out_pdf: str | Path, *, view: str = "both", pick: str = "both",
                 label_column: str | None = "Genes", max_hit_panels: int = 60) -> Path:
    """Fold-change and -log10(p) profiles across temperatures, per contrast.

    view: "FC" | "pval" | "both";  pick: "All_prot" | "Over_th" | "both" (R names).
    """
    if view not in ("FC", "pval", "both") or pick not in ("All_prot", "Over_th", "both"):
        raise ValueError("view must be FC/pval/both and pick All_prot/Over_th/both")
    temps = res.temperatures
    fc, lp = res.call_params.fc_cutoff, res.call_params.logp_cutoff
    metrics = []
    if view in ("FC", "both"):
        metrics.append(("diff", "log2 fold change", [0.0, fc, -fc]))
    if view in ("pval", "both"):
        metrics.append(("logpval", f"-log10({res.call_params.p_column})", [0.0, lp]))
    out_pdf = Path(out_pdf)
    labels_all = _label_series(res.data, label_column)
    lab_map = dict(zip(res.data["name"], labels_all))

    with _style(), PdfPages(out_pdf) as pdf:
        for c in res.contrasts:
            h = res.hits[res.hits["contrast"] == c]
            dir_map = dict(zip(h["name"], h["direction"]))
            groups = {n: g for n, g in res.data.groupby("name", sort=False)}

            if pick in ("All_prot", "both"):
                fig, axes = plt.subplots(1, len(metrics), figsize=(5.2 * len(metrics), 4.2), squeeze=False)
                for ax, (suffix, ylab, lines) in zip(axes[0], metrics):
                    col = f"{c}_{suffix}"
                    # background: one LineCollection (fast); a missing temperature leaves a gap
                    wide = res.data.pivot(index="name", columns="Temperature", values=col)
                    wide = wide.reindex(columns=list(temps.ids))
                    wide = wide[~wide.index.isin(list(dir_map))]
                    xs = np.array([temps.x(t) for t in temps.ids])
                    segs = [np.column_stack([xs, row]) for row in wide.to_numpy(float)]
                    ax.add_collection(LineCollection(segs, colors=NEUTRAL, linewidths=0.4,
                                                     alpha=0.35, zorder=1, rasterized=True))
                    ax.autoscale_view()
                    for n, d in dir_map.items():
                        x, y = _series_arrays(groups[n], col, temps)
                        color = UP if d == "up" else DOWN if d == "down" else INK_2
                        ax.plot(x, y, color=color, lw=1.2, marker="o", ms=2.5, zorder=3)
                        if len(h) <= 12 and np.isfinite(y).any():
                            i = int(np.nanargmax(np.abs(y)))
                            ax.annotate(lab_map.get(n, n), (x[i], y[i]), fontsize=6, color=INK,
                                        xytext=(3, 3), textcoords="offset points")
                    for v in lines:
                        ax.axhline(v, color=INK_2 if v == 0 else MUTED, lw=0.7,
                                   ls="-" if v == 0 else "--", zorder=2.5)
                    _set_temperature_axis(ax, temps)
                    ax.set_ylabel(ylab)
                fig.suptitle(f"{c} — all retained proteins ({len(groups)}), {len(h)} hits highlighted",
                             fontsize=10, fontweight="bold")
                fig.legend(handles=[Line2D([], [], color=UP, lw=1.2, label="Hit, up"),
                                    Line2D([], [], color=DOWN, lw=1.2, label="Hit, down"),
                                    Line2D([], [], color=NEUTRAL, lw=1, label="Other proteins")],
                           loc="upper center", ncol=3, fontsize=7, bbox_to_anchor=(0.5, 0.94))
                fig.tight_layout(rect=(0, 0, 1, 0.88))
                pdf.savefig(fig, dpi=200)
                plt.close(fig)

            if pick in ("Over_th", "both") and len(h):
                names = h.sort_values("max_logp", ascending=False)["name"].tolist()[:max_hit_panels]
                for suffix, ylab, lines in metrics:
                    col = f"{c}_{suffix}"
                    ncol = min(4, len(names))
                    per_page = ncol * 4
                    for start in range(0, len(names), per_page):
                        chunk = names[start:start + per_page]
                        nrow = math.ceil(len(chunk) / ncol)
                        fig, axes = plt.subplots(nrow, ncol, figsize=(2.7 * ncol + 0.4, 2.2 * nrow + 0.8),
                                                 sharex=True, sharey=True, squeeze=False)
                        for ax, n in zip(axes.flat, chunk):
                            x, y = _series_arrays(groups[n], col, temps)
                            _, dd = _series_arrays(groups[n], f"{c}_diff", temps)
                            _, ll = _series_arrays(groups[n], f"{c}_logpval", temps)
                            with np.errstate(invalid="ignore"):
                                passing = (np.abs(dd) >= fc) & (ll >= lp)
                            d = dir_map.get(n)
                            color = UP if d == "up" else DOWN if d == "down" else INK_2
                            ax.plot(x, y, color=color, lw=1.2, zorder=3)
                            ax.scatter(x[~passing], y[~passing], s=12, facecolors="white",
                                       edgecolors=color, linewidths=0.9, zorder=4)
                            ax.scatter(x[passing], y[passing], s=16, color=color, zorder=4)
                            for v in lines:
                                ax.axhline(v, color=INK_2 if v == 0 else MUTED, lw=0.6,
                                           ls="-" if v == 0 else "--")
                            row = h[h["name"] == n].iloc[0]
                            ax.set_title(f"{lab_map.get(n, n)} ({n})", fontsize=7)
                            ax.text(0.02, 0.96, f"{row['run_start']}–{row['run_end']} · {row['shift_class']}",
                                    transform=ax.transAxes, fontsize=5.5, color=INK_2, va="top")
                        for ax in list(axes.flat)[len(chunk):]:
                            ax.set_visible(False)
                        for ax in axes[-1, :]:
                            _set_temperature_axis(ax, temps)
                        for ax in axes[:, 0]:
                            ax.set_ylabel(ylab)
                        fig.suptitle(f"{c} — hits ({ylab}); filled = temperature passing both thresholds",
                                     fontsize=9, fontweight="bold")
                        fig.tight_layout(rect=(0, 0, 1, 0.94))
                        pdf.savefig(fig)
                        plt.close(fig)
    log.info("Profiles: %s", out_pdf)
    return out_pdf


def shift_pdf(sr, res: TPPResult, out_pdf: str | Path, *, label_column: str | None = "Genes",
              top: int = 12) -> Path:
    """Shift test figures, per contrast: p-value histogram (calibration check) and
    the top candidates' Δ profile (± 95 % CI) with baseline and best window."""
    temps = res.temperatures
    xs = np.array([temps.x(t) for t in temps.ids])
    out_pdf = Path(out_pdf)
    lab_map = dict(zip(res.data["name"], _label_series(res.data, label_column)))
    with _style(), PdfPages(out_pdf) as pdf:
        for j, c in enumerate(res.contrasts):
            tab = sr.table[sr.table["contrast"] == c].sort_values("p_protein")
            fig, ax = plt.subplots(figsize=(5.2, 3.6))
            pv = tab["p_protein"].dropna()
            ax.hist(pv, bins=np.linspace(0, 1, 21), color=NEUTRAL, edgecolor="white", linewidth=0.6)
            ax.axhline(len(pv) / 20, color=INK_2, lw=0.7, ls="--")
            ax.set_xlabel("protein-level p-value (Simes over windows)")
            ax.set_ylabel("proteins")
            ax.set_title(f"{c} — {len(pv)} proteins tested, {int(tab['hit'].sum())} hits "
                         f"(q ≤ {sr.params.alpha:g})", fontsize=9)
            ax.text(0.98, 0.95, "dashed: uniform (no signal)", transform=ax.transAxes, ha="right",
                    va="top", fontsize=6.5, color=INK_2)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

            sel = tab.head(top)
            if sel.empty:
                continue
            ncol = 4
            nrow = int(np.ceil(len(sel) / ncol))
            fig, axes = plt.subplots(nrow, ncol, figsize=(2.7 * ncol + 0.4, 2.3 * nrow + 0.9),
                                     sharex=True, squeeze=False)
            for ax, (_, row) in zip(axes.flat, sel.iterrows()):
                n = row["name"]
                d, se = [], []
                for t in temps.ids:
                    f = sr.fits[t]
                    k = f.proteins.get_indexer([n])[0]
                    d.append(f.diff[k, j] if k >= 0 else np.nan)
                    se.append(f.se[k, j] if k >= 0 else np.nan)
                d, se = np.array(d), np.array(se)
                color = UP if row["direction"] == "up" else DOWN if row["direction"] == "down" else INK_2
                if row["best_window"]:
                    parts = str(row["best_window"]).split("-")
                    x0, x1 = temps.x(parts[0]), temps.x(parts[-1])
                    pad = (xs[1] - xs[0]) / 2 if len(xs) > 1 else 0.5
                    ax.axvspan(x0 - pad, x1 + pad, color=color, alpha=0.10, lw=0)
                ax.axhline(0, color=INK_2, lw=0.6)
                if np.isfinite(row["baseline"]):
                    ax.axhline(row["baseline"], color=MUTED, lw=0.7, ls="--")
                ax.errorbar(xs, d, yerr=1.96 * se, color=color, lw=1.1, marker="o", ms=2.8,
                            elinewidth=0.6, capsize=1.5)
                star = " ★" if row["hit"] else ""
                ax.set_title(f"{lab_map.get(n, n)} ({n}){star}", fontsize=7)
                ax.text(0.02, 0.96, f"q={row['q_protein']:.2g} · {row['best_window']} · D={row['dd_best']:+.2f}",
                        transform=ax.transAxes, fontsize=5.5, color=INK_2, va="top")
            for ax in list(axes.flat)[len(sel):]:
                ax.set_visible(False)
            for ax in axes[-1, :]:
                _set_temperature_axis(ax, temps)
            for ax in axes[:, 0]:
                ax.set_ylabel("log2 FC (± 95 % CI)")
            fig.suptitle(f"{c} — top {len(sel)} shift-test candidates (dashed: baseline, shaded: best window, "
                         f"★ hit)", fontsize=9, fontweight="bold")
            fig.tight_layout(rect=(0, 0, 1, 0.94))
            pdf.savefig(fig)
            plt.close(fig)
    log.info("Shift test plots: %s", out_pdf)
    return out_pdf
