"""Command line: ``fctpp <command>`` (or ``python -m fctpp <command>``).

  fctpp init-config [-o fctpp_config.yaml]
  fctpp prepare     --raw-dir TPP_data --design TPP_data/experimental_design.xlsx --runs deimos_runs
  fctpp run-deimos  --runs deimos_runs --deimos-dir ../deimos-proteomics-pipeline [--jobs 4]
  fctpp analyze     --deimos-root deimos_runs --out fctpp_output
  fctpp all         (prepare + run-deimos + analyze)

Every option can live in the YAML (``--config``); CLI values win.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .config import load_config, parse_override, set_dotted

log = logging.getLogger("fctpp")


def _common(p: argparse.ArgumentParser) -> None:
    p.add_argument("-c", "--config", help="YAML configuration file")
    p.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                   help="override any key, e.g. --set call.fc_cutoff=0.8 (repeatable)")
    p.add_argument("--celsius", help="temperature → °C mapping (CSV: temperature,celsius)")
    p.add_argument("-v", "--verbose", action="store_true")


def _analysis_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--deimos-root", help="folder with one Deimos output per temperature")
    p.add_argument("--manifest", help="runs.csv (temperature, results_path[, temperature_c])")
    p.add_argument("--long-table", help="long table instead of Deimos workbooks (CSV/XLSX)")
    p.add_argument("--sheet", help="sheet of --long-table")
    p.add_argument("-o", "--out", help="output folder")
    p.add_argument("--fc", type=float, help="|log2 FC| threshold (default 0.6)")
    p.add_argument("--logp", type=float, help="-log10(p) threshold (default 2)")
    p.add_argument("--rule", choices=["strict", "r_legacy"], help="call rule (default strict)")
    p.add_argument("--min-temperatures", type=int)
    p.add_argument("--min-streak", type=int)
    p.add_argument("--no-plots", action="store_true", help="tables only")


def _prepare_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--raw-dir", help="folder with one DIA-NN pg_matrix per temperature (T1.tsv...)")
    p.add_argument("--design", help="experimental design (xlsx/csv)")
    p.add_argument("--runs", help="output folder for the Deimos runs (default deimos_runs)")


def _deimos_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--deimos-dir", help="path to the deimos-proteomics-pipeline repository")
    p.add_argument("--jobs", type=int, help="Deimos runs in parallel")
    p.add_argument("--skip-existing", action="store_true", help="do not re-run finished temperatures")
    if not any(a.dest == "runs" for a in p._actions):
        p.add_argument("--runs", help="folder written by 'fctpp prepare'")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="fctpp", description="FC.TPP — fold-change Thermal Proteome "
                                 "Profiling on per-temperature Deimos results.")
    ap.add_argument("--version", action="version", version=f"fctpp {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("init-config", help="write an annotated YAML template")
    p.add_argument("-o", "--output", default="fctpp_config.yaml")

    p = sub.add_parser("prepare", help="write one Deimos input set per temperature")
    _common(p)
    _prepare_args(p)

    p = sub.add_parser("run-deimos", help="run Deimos on every prepared temperature")
    _common(p)
    _deimos_args(p)

    p = sub.add_parser("analyze", help="TPP calls, Excel and figures from Deimos outputs")
    _common(p)
    _analysis_args(p)

    p = sub.add_parser("all", help="prepare + run-deimos + analyze")
    _common(p)
    _prepare_args(p)
    _deimos_args(p)
    _analysis_args(p)
    return ap


def _apply_cli(cfg: dict, a: argparse.Namespace) -> dict:
    def put(attr, key, transform=lambda v: v):
        v = getattr(a, attr, None)
        if v not in (None, False):
            set_dotted(cfg, key, transform(v))

    put("celsius", "temperatures_celsius")
    put("raw_dir", "prepare.raw_dir")
    put("design", "prepare.design")
    put("runs", "prepare.out_dir")
    put("deimos_dir", "prepare.deimos_dir")
    put("jobs", "prepare.jobs")
    put("skip_existing", "prepare.skip_existing")
    put("out", "output.out_dir")
    put("fc", "call.fc_cutoff")
    put("logp", "call.logp_cutoff")
    put("rule", "call.rule")
    put("min_temperatures", "filter.min_temperatures")
    put("min_streak", "filter.min_streak")
    if getattr(a, "manifest", None):
        set_dotted(cfg, "input.mode", "deimos")
        set_dotted(cfg, "input.manifest", a.manifest)
    if getattr(a, "deimos_root", None):
        set_dotted(cfg, "input.mode", "deimos")
        set_dotted(cfg, "input.deimos_root", a.deimos_root)
        set_dotted(cfg, "input.manifest", None)
    if getattr(a, "long_table", None):
        set_dotted(cfg, "input.mode", "long-table")
        set_dotted(cfg, "input.long_table", a.long_table)
        if a.sheet:
            set_dotted(cfg, "input.long_table_sheet", a.sheet)
    if getattr(a, "no_plots", False):
        for k in ("volcano", "profiles", "plot3d", "gif"):
            set_dotted(cfg, f"plots.{k}", False)
    for item in getattr(a, "set", []) or []:
        k, v = parse_override(item)
        set_dotted(cfg, k, v)
    return cfg


def main(argv: list[str] | None = None) -> int:
    a = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if getattr(a, "verbose", False) else logging.INFO,
                        format="%(levelname)-7s %(message)s")
    for noisy in ("fontTools", "matplotlib", "PIL"):
        logging.getLogger(noisy).setLevel(logging.WARNING)
    if a.command == "init-config":
        from .config import example_yaml
        Path(a.output).write_text(example_yaml(), encoding="utf-8")
        print(f"Template written: {a.output}")
        return 0

    cfg = _apply_cli(load_config(a.config), a)
    from . import pipeline

    try:
        if a.command in ("prepare", "all"):
            pipeline.prepare(cfg)
        if a.command in ("run-deimos", "all"):
            pipeline.run_deimos(cfg)
        if a.command == "all" and not cfg["input"].get("manifest") and not getattr(a, "deimos_root", None):
            cfg["input"]["mode"] = "deimos"
            cfg["input"]["manifest"] = str(Path(cfg["prepare"]["out_dir"]) / "runs.csv")
        if a.command in ("analyze", "all"):
            out = pipeline.analyze(cfg)
            res = out["result"]
            print("\nSummary")
            print(res.summary.to_string(index=False))
            print(f"\nResults in {Path(cfg['output']['out_dir']).resolve()}")
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        log.error("%s", e)
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
