"""High-level steps used by the CLI (and usable from Python)."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .config import save_config
from .export import write_workbook
from .inputs import (DEFAULT_WORKBOOK_GLOB, discover_deimos_runs, load_deimos_runs,
                     load_long_table, read_manifest)
from .temperatures import TemperatureSeries, read_celsius_mapping
from .tpp import CallParams, FilterParams, TPPResult, run_tpp

log = logging.getLogger(__name__)


def _safe(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in name)


def load_inputs(cfg: dict):
    i = cfg["input"]
    celsius = read_celsius_mapping(cfg.get("temperatures_celsius"))
    mode = i.get("mode", "deimos")
    if mode == "deimos":
        if i.get("manifest"):
            runs, cel2 = read_manifest(i["manifest"])
            celsius = celsius or cel2
            source = str(i["manifest"])
        else:
            runs = discover_deimos_runs(i["deimos_root"], i.get("temperature_regex"),
                                        i.get("workbook_glob") or DEFAULT_WORKBOOK_GLOB)
            source = str(i["deimos_root"])
        loaded = load_deimos_runs(runs)
    elif mode == "long-table":
        if not i.get("long_table"):
            raise ValueError("input.long_table is required in long-table mode")
        loaded = load_long_table(i["long_table"], i.get("long_table_sheet"))
        source = str(i["long_table"])
    else:
        raise ValueError(f"Unknown input.mode {mode!r} (deimos | long-table)")
    temps = TemperatureSeries.from_ids(loaded.results["Temperature"].unique(), celsius)
    if temps.celsius:
        loaded.results["Temperature_C"] = loaded.results["Temperature"].map(temps.celsius)
    return loaded, temps, mode, source


def analyze(cfg: dict) -> dict:
    """Load inputs → TPP calls → Excel/CSV → figures. Returns output paths."""
    out_dir = Path(cfg["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    loaded, temps, mode, source = load_inputs(cfg)
    log.info("%d temperatures: %s", len(temps), ", ".join(temps.ids))
    log.info("Contrasts: %s", ", ".join(loaded.contrasts))

    f, c = cfg["filter"], cfg["call"]
    res: TPPResult = run_tpp(
        loaded.results, temps,
        FilterParams(float(f["min_temperatures"]), float(f["min_streak"])),
        CallParams(fc_cutoff=float(c["fc_cutoff"]), logp_cutoff=float(c["logp_cutoff"]),
                   p_column=str(c["p_column"]), min_consecutive=int(c["min_consecutive"]),
                   rule=str(c["rule"]), abundance_check=bool(c.get("abundance_check", True))),
        contrasts=loaded.contrasts)
    res.notes.extend(loaded.warnings)

    outputs: dict[str, object] = {}
    from .shift_pipeline import annotate_legacy_hits, methods_rows, run_shift_analysis, summary_columns
    sr = run_shift_analysis(cfg, loaded, temps, res, mode)
    outputs["shift"] = sr
    extra, mrows = [], []
    if sr is not None:
        extra = [("Shift_test", sr.table)]
        if len(sr.fdr):
            extra.append(("FDR_estimates", sr.fdr))
        if sr.shared_table is not None and len(sr.shared_table):
            extra.append(("Shared_component", sr.shared_table))
        mrows = methods_rows(sr)
        res.hits = annotate_legacy_hits(sr, res.hits)
    run_info = {"input_mode": mode, "input": source, "warnings": loaded.warnings, "config": cfg}
    outputs["excel"] = write_workbook(res, out_dir / cfg["output"]["excel"], run_info=run_info,
                                      provenance=loaded.provenance, methods_upstream=loaded.methods,
                                      summary=summary_columns(sr, res), extra_sheets=extra,
                                      methods_extra=mrows)
    if cfg["output"].get("csv", True):
        res.data.to_csv(out_dir / "TPP_data.csv", index=False)
        res.hits.to_csv(out_dir / "TPP_hits.csv", index=False)
        outputs["csv"] = [out_dir / "TPP_data.csv", out_dir / "TPP_hits.csv"]
        if sr is not None:
            sr.table.to_csv(out_dir / "shift_test.csv", index=False)
            sr.windows.to_csv(out_dir / "shift_windows.csv", index=False)
            outputs["csv"] += [out_dir / "shift_test.csv", out_dir / "shift_windows.csv"]

    pl = cfg["plots"]
    label = pl.get("label_column")
    if pl.get("volcano", True):
        from .plots import volcano_pdf
        src = None
        if pl.get("volcano_all_proteins", True):
            src = loaded.results.copy()
            for ct in res.contrasts:
                p = pd.to_numeric(src[f"{ct}_{res.call_params.p_column}"], errors="coerce").astype(float)
                with np.errstate(divide="ignore"):
                    src[f"{ct}_logpval"] = -np.log10(p)
        outputs["volcano"] = volcano_pdf(res, out_dir / "volcano_plots.pdf",
                                         fc_cutoff=pl.get("volcano_fc_cutoff"),
                                         logp_cutoff=pl.get("volcano_logp_cutoff"),
                                         layout=pl.get("volcano_layout", "grid"),
                                         label_column=label, source=src)
    if pl.get("profiles", True):
        from .plots import profiles_pdf
        outputs["profiles"] = profiles_pdf(res, out_dir / "TPP_plots.pdf",
                                           view=pl.get("profiles_view", "both"),
                                           pick=pl.get("profiles_pick", "both"), label_column=label)
    if sr is not None and pl.get("profiles", True):
        from .plots import shift_pdf
        outputs["shift_plots"] = shift_pdf(sr, res, out_dir / "shift_test_plots.pdf", label_column=label)
    if pl.get("plot3d", True) or pl.get("gif", False):
        from .plot3d import plot3d_html, rotation_gif
        outputs["plot3d"], outputs["gif"] = [], []
        for ct in res.contrasts:
            if pl.get("plot3d", True):
                outputs["plot3d"].append(plot3d_html(res, ct, out_dir / f"{_safe(ct)}_plot3D.html",
                                                     label_column=label,
                                                     include_plotlyjs=pl.get("plotly_js", "embed")))
            if pl.get("gif", False):
                outputs["gif"].append(rotation_gif(res, ct, out_dir / f"{_safe(ct)}_rotation.gif",
                                                   label_column=label))
    save_config(cfg, out_dir / "last_config.yaml")
    outputs["result"] = res
    return outputs
