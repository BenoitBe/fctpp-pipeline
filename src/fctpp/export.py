"""Excel export — one workbook written once (no cross-function overwrite).

In R, ``format_tsv``, ``sumXP``, ``process_and_export`` and ``tidy_tpp`` all
defaulted to ``TPP_results.xlsx``; ``tidy_tpp`` recreated the file, so the
``comparaison_*`` sheets were lost (the delivered TPP_results.xlsx only holds
``TPP_data`` and ``Stats_proteins``).
"""

from __future__ import annotations

import json
import logging
import platform
import re
from datetime import datetime
from pathlib import Path

import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from . import __version__
from .tpp import TPPResult, wide_matrix

log = logging.getLogger(__name__)

HEADER_FILL = PatternFill("solid", fgColor="D9D9D9")   # platform convention (R/Deimos)
HEADER_FONT = Font(bold=True)


def _sheet_name(base: str, used: set[str]) -> str:
    name = re.sub(r"[\[\]:*?/\\]", "_", base)[:31]
    k, out = 1, name
    while out in used:
        suffix = f"~{k}"
        out = name[: 31 - len(suffix)] + suffix
        k += 1
    used.add(out)
    return out


def _style_sheet(ws, df: pd.DataFrame, freeze: str | None = "A2") -> None:
    for cell in ws[1]:
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(vertical="center")
    if freeze and ws.max_row > 1:
        ws.freeze_panes = freeze
    if ws.max_row > 1 and ws.max_column > 0:
        ws.auto_filter.ref = ws.dimensions
    for i, col in enumerate(df.columns, start=1):
        # pandas 3: astype(str) keeps NaN as NaN -> measure str() explicitly
        sample = df[col].head(200).map(lambda v: len(str(v)))
        width = max([len(str(col))] + sample.tolist()) if len(sample) else len(str(col))
        ws.column_dimensions[get_column_letter(i)].width = min(max(8, width + 2), 60)


def methods_table(res: TPPResult, run_info: dict) -> pd.DataFrame:
    cp, fp = res.call_params, res.filter_params
    rule = ("strict: >= {k} ADJACENT temperatures of the series, all with |log2FC| >= {fc} and "
            "-log10({p}) >= {lp}, with the SAME sign of fold change"
            if cp.rule == "strict" else
            "r_legacy (FC.TPP R): two SUCCESSIVE rows of the protein (temperature order, gaps "
            "allowed) with |log2FC| >= {fc} and -log10({p}) >= {lp}; sign not checked")
    rows = [
        ("Tool", f"FC.TPP (Python) {__version__}", ""),
        ("Run date", datetime.now().strftime("%Y-%m-%d %H:%M"), ""),
        ("Python", platform.python_version(), f"pandas {pd.__version__}"),
        ("Input mode", run_info.get("input_mode", ""), run_info.get("input", "")),
        ("Upstream statistics", "Deimos, one run per temperature",
         "Filtering, imputation and limma/eBayes per temperature — see the Inputs sheet "
         "(Deimos 'Methods' values per temperature)"),
        ("Temperatures", f"{len(res.temperatures)}", ", ".join(res.temperatures.label(t) for t in res.temperatures)),
        ("Contrasts", f"{len(res.contrasts)}", ", ".join(res.contrasts)),
        ("Presence filter", f"nb_temp >= {fp.min_temperatures} and max_consecutive >= {fp.min_streak}",
         "nb_temp = temperatures where the protein passed Deimos' filter; max_consecutive = "
         "longest run of adjacent temperatures"),
        ("Proteins retained", f"{int(res.stats['retained'].sum())} / {len(res.stats)}",
         f"{res.n_lost} proteins lost after filtering"),
        ("P-value used", cp.p_column, "-log10 stored in <contrast>_logpval"),
        ("Call rule", cp.rule, rule.format(k=cp.min_consecutive, fc=cp.fc_cutoff,
                                           p=cp.p_column, lp=cp.logp_cutoff)),
        ("Abundance check", str(cp.abundance_check),
         f"hit also passing the thresholds with the same sign at the reference temperature "
         f"({res.temperatures.reference}) → 'abundance-like' (difference present without heating; "
         f"expression/solubility rather than thermal stability)"),
        ("Imputation flag", "imputed_in_run",
         "True when Deimos imputed at least one value of the protein (any condition: Deimos "
         "num_NAs is per protein) at a temperature of the calling run: check the hit on the raw data"),
        ("Multiple testing", "not controlled across temperatures",
         "Per-temperature p-values are thresholded; the consecutive-temperature rule reduces "
         "but does not calibrate false positives (see README)"),
    ]
    for w in run_info.get("warnings", []):
        rows.append(("Warning", "", w))
    return pd.DataFrame(rows, columns=["Item", "Value", "Detail"])


def write_workbook(res: TPPResult, out_path: str | Path, *, run_info: dict | None = None,
                   provenance: pd.DataFrame | None = None, methods_upstream: pd.DataFrame | None = None,
                   wide: bool = True) -> Path:
    run_info = run_info or {}
    out_path = Path(out_path)
    used: set[str] = set()
    sheets: list[tuple[str, pd.DataFrame]] = [
        ("Methods", methods_table(res, run_info)),
        ("Summary", res.summary),
        ("Hits", res.hits),
        ("TPP_data", res.data),
        ("Stats_proteins", res.stats),
        ("Calls", res.calls),
    ]
    if wide:
        for c in res.contrasts:
            sheets.append((f"FC_{c}", wide_matrix(res, c, "diff").reset_index()))
            sheets.append((f"logP_{c}", wide_matrix(res, c, "logpval").reset_index()))
    if provenance is not None and len(provenance):
        sheets.append(("Inputs", provenance))
    if methods_upstream is not None and len(methods_upstream):
        sheets.append(("Deimos_Methods", methods_upstream))
    sheets.append(("Temperatures", res.temperatures.to_frame()))
    sheets.append(("Parameters", pd.DataFrame(
        [("fctpp_params", json.dumps(run_info.get("config", res.params_dict()), default=str, indent=1))],
        columns=["key", "value"])))

    with pd.ExcelWriter(out_path, engine="openpyxl") as xw:
        for base, df in sheets:
            name = _sheet_name(base, used)
            df.to_excel(xw, sheet_name=name, index=False)
            _style_sheet(xw.sheets[name], df, freeze="C2" if base in ("TPP_data",) else "A2")
        ws = xw.sheets["Methods"]
        ws.column_dimensions["C"].width = 110
        for row in ws.iter_rows(min_row=2):
            row[2].alignment = Alignment(wrap_text=True, vertical="top")
    log.info("Excel: %s", out_path)
    return out_path
