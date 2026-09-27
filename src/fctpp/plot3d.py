"""Interactive 3D view (plotly) and optional rotating GIF (matplotlib).

Port of the 3D part of ``plot_and_3D()``. Differences with R:
* non-hit proteins are packed in ONE trace (NaN-separated lines) instead of one
  trace per protein, which keeps the HTML small and fast;
* the threshold planes span the actual temperature series (R hard-coded T1/T12);
* the GIF is rendered with matplotlib (no webshot2/Chrome/ImageMagick needed).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .plots import DOWN, INK, INK_2, NEUTRAL, UP, _label_series
from .tpp import TPPResult

log = logging.getLogger(__name__)


def _packed_lines(df: pd.DataFrame, c: str, temps) -> tuple[list, list, list, list]:
    xs, ys, zs, texts = [], [], [], []
    for name, g in df.groupby("name", sort=False):
        g = g.sort_values("Temperature", key=lambda s: s.map(temps.position))
        xs += [temps.x(t) for t in g["Temperature"]] + [None]
        ys += pd.to_numeric(g[f"{c}_diff"], errors="coerce").round(3).tolist() + [None]
        zs += pd.to_numeric(g[f"{c}_logpval"], errors="coerce").round(3).tolist() + [None]
        texts += [name] * len(g) + [None]
    return xs, ys, zs, texts


def plot3d_html(res: TPPResult, contrast: str, out_html: str | Path, *,
                label_column: str | None = "Genes", include_plotlyjs: str = "embed") -> Path:
    import plotly.graph_objects as go

    temps = res.temperatures
    fc, lp = res.call_params.fc_cutoff, res.call_params.logp_cutoff
    c = contrast
    data = res.data
    h = res.hits[res.hits["contrast"] == c]
    dir_map = dict(zip(h["name"], h["direction"]))
    labels = dict(zip(data["name"], _label_series(data, label_column)))

    fig = go.Figure()
    under = data[~data["name"].isin(dir_map)]
    xs, ys, zs, texts = _packed_lines(under, c, temps)
    fig.add_trace(go.Scatter3d(
        x=xs, y=ys, z=zs, text=texts, mode="lines+markers", name=f"Under threshold ({under['name'].nunique()})",
        # WebGL 3D lines ignore rgba alpha (they render white): use an opaque light grey
        line=dict(color="#d9d8d2", width=1),
        marker=dict(color="#c3c2b7", size=1.5),
        hovertemplate="%{text}<br>x=%{x}<br>log2FC=%{y:.2f}<br>-log10p=%{z:.2f}<extra></extra>"))

    many = len(dir_map) > 40
    for d, color in (("up", UP), ("down", DOWN), ("mixed", INK_2)):
        names = [n for n, dd in dir_map.items() if dd == d]
        if not names:
            continue
        if many:
            sub = data[data["name"].isin(names)]
            xs, ys, zs, texts = _packed_lines(sub, c, temps)
            fig.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, text=texts, mode="lines+markers",
                                       name=f"Hits {d} ({len(names)})",
                                       line=dict(color=color, width=6), marker=dict(color=color, size=3),
                                       hovertemplate="%{text}<br>log2FC=%{y:.2f}<br>-log10p=%{z:.2f}<extra></extra>"))
            continue
        for n in names:
            g = data[data["name"] == n].sort_values("Temperature", key=lambda s: s.map(temps.position))
            fig.add_trace(go.Scatter3d(
                x=[temps.x(t) for t in g["Temperature"]], y=g[f"{c}_diff"], z=g[f"{c}_logpval"],
                mode="lines+markers", name=f"{labels.get(n, n)} ({d})", legendgroup=d,
                line=dict(color=color, width=6), marker=dict(color=color, size=4),
                hovertemplate=f"{labels.get(n, n)} ({n})<br>log2FC=%{{y:.2f}}<br>-log10p=%{{z:.2f}}<extra></extra>"))

    # label at the most significant temperature of each hit (R: slice_max(logpval))
    if dir_map:
        sub = data[data["name"].isin(dir_map)]
        top = sub.loc[sub.groupby("name")[f"{c}_logpval"].idxmax()]
        fig.add_trace(go.Scatter3d(
            x=[temps.x(t) for t in top["Temperature"]], y=top[f"{c}_diff"], z=top[f"{c}_logpval"],
            mode="text", text=[labels.get(n, n) for n in top["name"]], textposition="top center",
            textfont=dict(size=10, color=INK), showlegend=False, hoverinfo="skip"))

    x0, x1 = temps.x(temps.ids[0]), temps.x(temps.ids[-1])
    yv = pd.to_numeric(data[f"{c}_diff"], errors="coerce")
    zv = pd.to_numeric(data[f"{c}_logpval"], errors="coerce").replace([np.inf], np.nan)
    ymin, ymax = float(np.nanmin([yv.min(), -fc])), float(np.nanmax([yv.max(), fc]))
    zmax = float(np.nanmax([zv.max(), lp]))
    plane = dict(color="rgb(137,135,129)", opacity=0.10, showscale=False, hoverinfo="skip",
                 showlegend=False)
    for yy in (fc, -fc):   # vertical planes |log2FC| = fc
        fig.add_trace(go.Mesh3d(x=[x0, x1, x1, x0], y=[yy] * 4, z=[0, 0, zmax, zmax],
                                i=[0, 0], j=[1, 2], k=[2, 3], **plane))
    fig.add_trace(go.Mesh3d(x=[x0, x1, x1, x0], y=[ymin, ymin, ymax, ymax], z=[lp] * 4,
                            i=[0, 0], j=[1, 2], k=[2, 3], **plane))

    axis_style = dict(backgroundcolor="#fcfcfb", gridcolor="#e1e0d9", zerolinecolor="#898781",
                      showbackground=True)
    fig.update_layout(
        title=f"{c} — {len(dir_map)} hits (|log2FC| ≥ {fc:g}, -log10p ≥ {lp:g}, "
              f"{res.call_params.min_consecutive} consecutive temperatures)",
        template="plotly_white",
        scene=dict(
            xaxis=dict(title=temps.axis_title, tickvals=[temps.x(t) for t in temps.ids],
                       ticktext=temps.tick_labels(), **axis_style),
            yaxis=dict(title="log2 fold change", **axis_style),
            zaxis=dict(title=f"-log10({res.call_params.p_column})", **axis_style),
            camera=dict(eye=dict(x=1.6, y=-1.6, z=0.8))),
        legend=dict(itemsizing="constant", font=dict(size=10)),
        margin=dict(l=0, r=0, t=50, b=0), paper_bgcolor="#fcfcfb")
    out_html = Path(out_html)
    inc = {"embed": True, "cdn": "cdn", "directory": "directory"}.get(include_plotlyjs, True)
    fig.write_html(str(out_html), include_plotlyjs=inc, full_html=True)
    log.info("3D plot: %s", out_html)
    return out_html


def rotation_gif(res: TPPResult, contrast: str, out_gif: str | Path, *, step_deg: int = 10,
                 fps: int = 10, label_column: str | None = "Genes",
                 max_background: int = 600, seed: int = 0) -> Path:
    """Rotating 3D view rendered with matplotlib (replaces webshot2 + magick).

    Rendering cost grows with the number of lines (~0.5 s per frame for 2,000
    proteins), so non-hit proteins are randomly subsampled to ``max_background``
    (fixed seed); hits are always drawn. The HTML view keeps every protein.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    temps = res.temperatures
    c = contrast
    data = res.data
    h = res.hits[res.hits["contrast"] == c]
    dir_map = dict(zip(h["name"], h["direction"]))
    labels = dict(zip(data["name"], _label_series(data, label_column)))

    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    fig = plt.figure(figsize=(7, 5.5), dpi=80)
    ax = fig.add_subplot(111, projection="3d")
    background = []                      # one collection instead of one artist per protein
    non_hits = [n for n in data["name"].unique() if n not in dir_map]
    if len(non_hits) > max_background:
        rng = np.random.default_rng(seed)
        keep = set(rng.choice(non_hits, size=max_background, replace=False)) | set(dir_map)
        data = data[data["name"].isin(keep)]
    for n, g in data.groupby("name", sort=False):
        g = g.sort_values("Temperature", key=lambda s: s.map(temps.position))
        x = np.array([temps.x(t) for t in g["Temperature"]])
        y = pd.to_numeric(g[f"{c}_diff"], errors="coerce").to_numpy(float)
        z = pd.to_numeric(g[f"{c}_logpval"], errors="coerce").to_numpy(float)
        d = dir_map.get(n)
        if d is None:
            pts = np.column_stack([x, y, z])
            pts = pts[np.isfinite(pts).all(axis=1)]
            if len(pts) > 1:            # 2-point segments: homogeneous array for mplot3d
                background.append(np.stack([pts[:-1], pts[1:]], axis=1))
        else:
            ax.plot(x, y, z, color=UP if d == "up" else DOWN if d == "down" else INK_2, lw=1.4)
            i = int(np.nanargmax(z))
            ax.text(x[i], y[i], z[i], labels.get(n, n), fontsize=6, color=INK)
    if background:
        segs = np.concatenate(background)
        ax.add_collection3d(Line3DCollection(segs, colors=NEUTRAL, linewidths=0.3, alpha=0.25))
        allp = segs.reshape(-1, 3)
        ax.set_xlim(allp[:, 0].min(), allp[:, 0].max())
        ax.set_ylim(np.nanmin(allp[:, 1]), np.nanmax(allp[:, 1]))
        ax.set_zlim(0, np.nanmax(allp[:, 2]) * 1.05)
    ax.set_xlabel(temps.axis_title, fontsize=7)
    ax.set_ylabel("log2 FC", fontsize=7)
    ax.set_zlabel("-log10 p", fontsize=7)
    ax.set_xticks([temps.x(t) for t in temps.ids])
    ax.set_xticklabels(temps.tick_labels(), fontsize=5)
    shown = f" (background: {min(len(non_hits), max_background)} of {len(non_hits)} proteins)" \
        if len(non_hits) > max_background else ""
    ax.set_title(f"{c} — {len(dir_map)} hits{shown}", fontsize=8)

    angles = list(range(0, 360, step_deg))

    def _update(a):
        ax.view_init(elev=22, azim=a)
        return ()

    anim = FuncAnimation(fig, _update, frames=angles, blit=False)
    out_gif = Path(out_gif)
    anim.save(str(out_gif), writer=PillowWriter(fps=fps))
    plt.close(fig)
    log.info("GIF: %s", out_gif)
    return out_gif
