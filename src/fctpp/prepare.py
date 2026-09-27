"""Prepare one Deimos run per temperature, and optionally run them.

Replaces the DEP half of the R package (``process_tsv`` → ``format_tsv`` →
``sumXP`` → ``process_and_export``): per-temperature filtering, imputation and
limma testing are delegated to Deimos, so the platform keeps a single
statistical engine.

What this module fixes on the way
---------------------------------
* Deimos reads the first **6** columns of a pg_matrix as annotations
  (``load_data``: ``tsv.columns[:6]``). Older DIA-NN exports (e.g. the TPP data:
  Protein.Group, Protein.Names, Genes, First.Protein.Description) only have 4,
  so the first two runs would be swallowed as annotations. Matrices are rewritten
  in the DIA-NN 2.x layout (6 annotation columns, missing ones left empty).
* Sample labels change with the temperature (``Sple-1B`` at T1, ``Sple-2B`` at
  T2...). A per-temperature ``ExperimentalDesign.csv`` is written with the exact
  labels Deimos will derive from the run names, and the mapping is traced.
* The design order is kept identical at every temperature, and the control
  condition can be forced, so that the contrast names/orientation are the same
  in the 12 Deimos outputs.
"""

from __future__ import annotations

import logging
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .design import DEIMOS_SAMPLE_REGEX, match_samples, read_design
from .temperatures import natural_key

log = logging.getLogger(__name__)

DIANN_ANNOTATION_COLUMNS = [
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
    "First.Protein.Description", "Proteotypic", "N.Sequences", "N.Proteotypic.Sequences",
]
# Layout written for Deimos (DIA-NN >= 2.0 pg_matrix): exactly 6 annotation columns.
DEIMOS_ANNOTATION_LAYOUT = [
    "Protein.Group", "Protein.Names", "Genes", "First.Protein.Description",
    "N.Sequences", "N.Proteotypic.Sequences",
]

DEFAULT_DEIMOS_PARAMS = {
    # statistics (see Deimos config_example.yaml)
    "impute_method": "qrilc",
    "fdr_global": False,
    "use_deqms": False,
    "n_iter_robustness": 0,
    "control_condition": "__auto__",
    "seed": 42,                 # honoured by Deimos versions with a seeded imputation
    # thresholds used by Deimos' own figures only (FC.TPP applies its own)
    "volcano_use_padj": False,
    "volcano_p_thresh": 0.01,
    "volcano_ratio_min": 1.5157,  # 2**0.6
    "anova_use_padj": False,
    "anova_p_thresh": 0.05,
    "n_heatmap_clusters": 3,
    # heavy / interactive modules off: 12 runs, no prompt allowed
    "go_organism": None,
    "run_perseverance": False,
    "run_gsea": False,
    "make_wgcna": False,
    "make_dashboard": False,
    "species_split": False,
}

RESULTS_WORKBOOK = "ProteomicAnalysis_Results.xlsx"


def discover_raw_files(folder: str | Path, pattern: str = "*.tsv") -> dict[str, Path]:
    """Temperature id (file stem) → path, naturally sorted (T1, T2, ..., T12)."""
    folder = Path(folder)
    if not folder.is_dir():
        raise FileNotFoundError(f"Raw data folder not found: {folder}")
    files = sorted(folder.glob(pattern), key=lambda p: natural_key(p.stem))
    if not files:
        raise FileNotFoundError(f"No '{pattern}' file in {folder}")
    return {p.stem: p for p in files}


def split_pg_matrix(df: pd.DataFrame) -> tuple[list[str], list[str]]:
    """(annotation columns, sample columns) of a DIA-NN pg_matrix."""
    ann = [c for c in df.columns if c in DIANN_ANNOTATION_COLUMNS]
    samples = [c for c in df.columns if c not in DIANN_ANNOTATION_COLUMNS]
    bad = [c for c in samples if not pd.api.types.is_numeric_dtype(df[c])]
    if bad:
        raise ValueError(f"Non-numeric sample column(s) {bad[:5]}: unknown annotation column?")
    if "Protein.Group" not in ann:
        raise ValueError("pg_matrix without a 'Protein.Group' column")
    return ann, samples


def to_deimos_layout(df: pd.DataFrame, sample_columns: list[str]) -> pd.DataFrame:
    """Rewrite a pg_matrix with the 6 annotation columns Deimos expects."""
    out = pd.DataFrame(index=df.index)
    for c in DEIMOS_ANNOTATION_LAYOUT:
        out[c] = df[c] if c in df.columns else np.nan
    for c in sample_columns:
        out[c] = df[c]
    return out


def _design_for_temperature(design: pd.DataFrame, temp: str) -> pd.DataFrame:
    if "temperature" in design.columns and design["temperature"].notna().any():
        sub = design[design["temperature"].astype(str) == str(temp)]
        if sub.empty:
            raise ValueError(f"No design row for temperature '{temp}' "
                             f"(temperatures in design: {sorted(design['temperature'].dropna().unique())})")
        return sub.reset_index(drop=True)
    return design


def prepare_deimos_inputs(raw_dir: str | Path, design_path: str | Path, out_dir: str | Path, *,
                          pattern: str = "*.tsv",
                          sample_regex: str | None = DEIMOS_SAMPLE_REGEX,
                          match_strategy: str = "auto",
                          deimos_params: dict | None = None,
                          celsius: dict[str, float] | None = None) -> pd.DataFrame:
    """Write ``<out_dir>/<T>/{report.pg_matrix.tsv, ExperimentalDesign.csv,
    deimos_config.yaml}`` for every raw file, plus ``runs.csv`` (manifest) and
    ``sample_mapping.csv``. Returns the manifest."""
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    raw = discover_raw_files(raw_dir, pattern)
    design = read_design(design_path)
    params = {**DEFAULT_DEIMOS_PARAMS, **(deimos_params or {})}

    manifest, mappings = [], []
    for temp, path in raw.items():
        df = pd.read_csv(path, sep="\t", low_memory=False)
        _, samples = split_pg_matrix(df)
        d_t = _design_for_temperature(design, temp)
        mapping = match_samples(samples, d_t, sample_regex=sample_regex, strategy=match_strategy)
        mapping.insert(0, "temperature", temp)
        mappings.append(mapping)

        unused = [c for c in samples if c not in set(mapping["column"])]
        if unused:
            log.warning("[%s] %d run column(s) not in the design, dropped: %s", temp, len(unused), unused)

        t_dir = out_dir / temp
        t_dir.mkdir(exist_ok=True)
        tsv_out = t_dir / "report.pg_matrix.tsv"
        to_deimos_layout(df, mapping["column"].tolist()).to_csv(tsv_out, sep="\t", index=False)

        des_out = t_dir / "ExperimentalDesign.csv"
        extra = [c for c in d_t.columns if c not in ("label", "condition", "replicate", "temperature")]
        des = pd.DataFrame({"label": mapping["sample_label"], "condition": mapping["condition"],
                            "replicate": mapping["replicate"]})
        for c in extra:          # e.g. subject / batch columns are passed through
            des[c] = d_t[c].to_numpy()
        des.to_csv(des_out, sep=";", index=False)

        cfg = dict(params)
        cfg.update({"tsv_path": str(tsv_out), "design_path": str(des_out),
                    "pr_path": str(t_dir / "report.pr_matrix.tsv"),
                    "out_dir": str(t_dir / "deimos_output")})
        cfg_out = t_dir / "deimos_config.yaml"
        with open(cfg_out, "w", encoding="utf-8") as fh:
            fh.write(f"# Deimos configuration generated by FC.TPP for temperature {temp}\n")
            yaml.safe_dump(cfg, fh, sort_keys=False, allow_unicode=True)

        manifest.append({
            "temperature": temp,
            "temperature_c": (celsius or {}).get(temp),
            "raw_file": str(Path(path).resolve()),
            "n_proteins_raw": len(df),
            "n_samples": len(mapping),
            "match": mapping["match"].iloc[0],
            "deimos_config": str(cfg_out),
            "results_path": str(t_dir / "deimos_output" / RESULTS_WORKBOOK),
        })
        log.info("[%s] %d proteins, %d samples (%s)", temp, len(df), len(mapping), mapping["match"].iloc[0])

    man = pd.DataFrame(manifest)
    man.to_csv(out_dir / "runs.csv", index=False)
    pd.concat(mappings, ignore_index=True).to_csv(out_dir / "sample_mapping.csv", index=False)
    log.info("Deimos inputs written for %d temperatures in %s", len(man), out_dir)
    return man


def _deimos_script(deimos_dir: str | Path) -> Path:
    p = Path(deimos_dir).expanduser().resolve()
    script = p / "deimos.py" if p.is_dir() else p
    if not script.exists():
        raise FileNotFoundError(f"deimos.py not found in {deimos_dir}")
    return script


def run_deimos_batch(manifest: pd.DataFrame | str | Path, deimos_dir: str | Path, *,
                     python: str | None = None, jobs: int = 1,
                     timeout: float | None = None, skip_existing: bool = False) -> pd.DataFrame:
    """Run ``python deimos.py --config <T>/deimos_config.yaml`` for each temperature.

    stdin is closed: if a Deimos version still asks a question, the run fails
    immediately (EOFError in its log) instead of hanging.
    """
    if not isinstance(manifest, pd.DataFrame):
        man_path = Path(manifest)
        if man_path.is_dir():
            man_path = man_path / "runs.csv"
        manifest = pd.read_csv(man_path)
    script = _deimos_script(deimos_dir)
    python = python or sys.executable

    def _one(row) -> dict:
        cfg = Path(row["deimos_config"])
        res = Path(row["results_path"])
        logf = cfg.parent / "deimos_log.txt"
        if skip_existing and res.exists():
            return {"temperature": row["temperature"], "status": "skipped (exists)", "seconds": 0.0,
                    "log": str(logf)}
        t0 = time.time()
        with open(logf, "w", encoding="utf-8") as fh:
            try:
                proc = subprocess.run([python, str(script), "--config", str(cfg)],
                                      cwd=str(cfg.parent), stdin=subprocess.DEVNULL,
                                      stdout=fh, stderr=subprocess.STDOUT, timeout=timeout,
                                      env={**_env(), "PYTHONIOENCODING": "utf-8", "MPLBACKEND": "Agg"})
                code = proc.returncode
            except subprocess.TimeoutExpired:
                code = "timeout"
        ok = code == 0 and res.exists()
        status = "ok" if ok else f"failed (exit {code})"
        log.info("[%s] Deimos %s in %.0f s", row["temperature"], status, time.time() - t0)
        return {"temperature": row["temperature"], "status": status,
                "seconds": round(time.time() - t0, 1), "log": str(logf)}

    rows = manifest.to_dict("records")
    if jobs and jobs > 1:
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            status = list(ex.map(_one, rows))
    else:
        status = [_one(r) for r in rows]
    st = pd.DataFrame(status)
    out = manifest.drop(columns=[c for c in ("status", "seconds", "log") if c in manifest.columns])
    out = out.merge(st, on="temperature", how="left")
    failed = out[~out["status"].astype(str).str.startswith(("ok", "skipped"))]
    if len(failed):
        log.error("Deimos failed for %s — see the deimos_log.txt files", failed["temperature"].tolist())
    return out


def _env() -> dict:
    import os
    return dict(os.environ)
