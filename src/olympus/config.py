"""YAML configuration (nested sections), with defaults and CLI overrides."""

from __future__ import annotations

import copy
import logging
from pathlib import Path
from typing import Any

import yaml

from .inputs import DEFAULT_WORKBOOK_GLOB

log = logging.getLogger(__name__)

DEFAULTS: dict[str, Any] = {
    "input": {
        "mode": "deimos",            # deimos | long-table
        "deimos_root": "deimos_results",  # folder holding the X Deimos result workbooks
        "workbook_glob": DEFAULT_WORKBOOK_GLOB,
        "temperature_regex": None,   # e.g. '(T\\d+)': temperature id taken from the path
        "manifest": None,            # CSV: temperature, results_path[, temperature_c]
        "long_table": None,          # CSV/XLSX in long format (e.g. R TPP_data sheet)
        "long_table_sheet": None,
    },
    # {T1: 37.0, T2: 40.4, ...} or a CSV path (temperature, celsius); null = ids only
    "temperatures_celsius": None,
    "filter": {                      # integer = number of temperatures,
        "min_temperatures": 6,       # value < 1 = fraction of the X temperatures
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
    "shift_test": {
        "enabled": "auto",           # auto = when Deimos matrices/ sit next to every workbook
        "deimos_dir": None,          # Deimos folder providing limma_ebayes.py (or env DEIMOS_DIR)
        "paired": "auto",            # auto (measured within-lysate correlation) | true | false
        "baseline": "none",          # none | median | lowest (see shift.py)
        "localisation_p": 0.01,      # detected window must differ from the profile median (null = off)
        "windows": [2, 3, 4],        # window sizes (adjacent temperatures), combined by Simes
        "window": None,              # single size instead (overrides windows); both null = call.min_consecutive
        "alpha": 0.05,               # protein-level FDR (Simes over windows, then BH)
        "fc_cutoff": None,           # null = call.fc_cutoff
        "use_adjacent_correlation": True,
        "shared_component": True,    # diagnostic of effects common to many proteins (flags hits)
        "shared_r_threshold": 0.7,   # |r| with the shared pattern above which a protein is flagged
        "shared_n_null": 100,        # sign-flip draws for the significance of the component
        "shared_min_ratio": 1.5,     # flags only if the component explains >= 1.5x its null 95 % quantile
        "permutations": True,        # within-lysate swaps -> FDR estimate of each rule
        "max_patterns": 64,
        "seed": 42,
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
        "out_dir": "olympus_output",
        "excel": "Olympus_Results.xlsx",
        "csv": True,
    },
}


def deep_merge(base: dict, over: dict | None, path: str = "") -> dict:
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if k not in out:
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
