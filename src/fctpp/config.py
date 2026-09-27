"""YAML configuration (nested sections), with defaults and CLI overrides."""

from __future__ import annotations

import copy
import logging
from pathlib import Path
from typing import Any

import yaml

from .design import DEIMOS_SAMPLE_REGEX
from .prepare import DEFAULT_DEIMOS_PARAMS

log = logging.getLogger(__name__)

DEFAULTS: dict[str, Any] = {
    "input": {
        "mode": "deimos",            # deimos | long-table
        "deimos_root": "deimos_runs",  # folder with one sub-folder per temperature
        "manifest": None,            # runs.csv (temperature, results_path[, temperature_c])
        "long_table": None,          # CSV/XLSX in long format (e.g. R TPP_data sheet)
        "long_table_sheet": None,
        "temperature_regex": None,   # keep only sub-folders matching this regex
    },
    # {T1: 37.0, T2: 40.4, ...} or a CSV path (temperature, celsius); null = ids only
    "temperatures_celsius": None,
    "filter": {
        "min_temperatures": 6,
        "min_streak": 6,
    },
    "call": {
        "fc_cutoff": 0.6,
        "logp_cutoff": 2.0,
        "p_column": "p.val",         # p.val (R behaviour) | p.adj
        "min_consecutive": 2,
        "rule": "strict",            # strict | r_legacy
        "abundance_check": True,
    },
    "plots": {
        "volcano": True,
        "volcano_layout": "grid",    # grid (1 page / contrast) | pages (1 page / temperature x contrast)
        "volcano_all_proteins": True,  # all quantified proteins (R) or only retained ones
        "volcano_fc_cutoff": None,   # null = call.fc_cutoff
        "volcano_logp_cutoff": None,  # null = call.logp_cutoff
        "profiles": True,
        "profiles_view": "both",     # FC | pval | both
        "profiles_pick": "both",     # All_prot | Over_th | both
        "plot3d": True,
        "plotly_js": "embed",        # embed (offline, ~3.5 MB/file) | cdn | directory
        "gif": False,
        "label_column": "Genes",
    },
    "output": {
        "out_dir": "fctpp_output",
        "excel": "FCTPP_Results.xlsx",
        "csv": True,
    },
    "prepare": {
        "raw_dir": "TPP_data",
        "pattern": "*.tsv",
        "design": None,              # xlsx/csv; template (R style) or long (with 'temperature')
        "out_dir": "deimos_runs",
        "sample_regex": DEIMOS_SAMPLE_REGEX,
        "match_strategy": "auto",    # auto | exact | suffix
        "deimos_dir": None,          # clone of deimos-proteomics-pipeline
        "python": None,              # interpreter for Deimos (default: current)
        "jobs": 1,
        "timeout": None,
        "skip_existing": False,
        "deimos": copy.deepcopy(DEFAULT_DEIMOS_PARAMS),
    },
}


def deep_merge(base: dict, over: dict | None, path: str = "") -> dict:
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if k not in out and path != "prepare.deimos":   # Deimos keys pass through freely
            log.warning("Unknown configuration key '%s%s' (kept)", f"{path}." if path else "", k)
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v, f"{path}.{k}" if path else k)
        else:
            out[k] = v
    return out


def load_config(path: str | Path | None = None) -> dict:
    cfg = copy.deepcopy(DEFAULTS)
    if path:
        with open(path, encoding="utf-8") as fh:
            user = yaml.safe_load(fh) or {}
        cfg = deep_merge(cfg, user)
    return cfg


def set_dotted(cfg: dict, dotted: str, value: Any) -> None:
    keys = dotted.split(".")
    node = cfg
    for k in keys[:-1]:
        node = node.setdefault(k, {})
    node[keys[-1]] = value


def parse_override(item: str) -> tuple[str, Any]:
    """'call.fc_cutoff=0.8' -> ('call.fc_cutoff', 0.8) (value parsed as YAML)."""
    if "=" not in item:
        raise ValueError(f"--set expects key=value, got {item!r}")
    k, v = item.split("=", 1)
    return k.strip(), yaml.safe_load(v)


def save_config(cfg: dict, path: str | Path) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False, allow_unicode=True)


def example_yaml() -> str:
    return (Path(__file__).with_name("config_template.yaml")).read_text(encoding="utf-8")
