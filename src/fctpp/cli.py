"""Command line: ``fctpp <command>`` (or ``python -m fctpp <command>``).

  fctpp init-config [-o fctpp_config.yaml]
  fctpp analyze --deimos-root deimos_results -o fctpp_output
  fctpp analyze --manifest runs.csv -o fctpp_output
  fctpp analyze --long-table TPP_results.xlsx --sheet TPP_data --rule r_legacy

Input = the X Deimos result workbooks (one per temperature). Every option can
live in the YAML (``--config``); CLI values win.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .config import load_config, parse_override, set_dotted

log = logging.getLogger("fctpp")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="fctpp", description="FC.TPP — fold-change Thermal Proteome "
                                 "Profiling on the per-temperature Deimos result workbooks.")
    ap.add_argument("--version", action="version", version=f"fctpp {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("init-config", help="write an annotated YAML template")
    p.add_argument("-o", "--output", default="fctpp_config.yaml")

    p = sub.add_parser("analyze", help="TPP calls, Excel report and figures")
    p.add_argument("-c", "--config", help="YAML configuration file")
    src = p.add_argument_group("inputs (one of)")
    src.add_argument("--deimos-root", help="folder holding the X Deimos result workbooks")
    src.add_argument("--manifest", help="CSV: temperature, results_path[, temperature_c]")
    src.add_argument("--long-table", help="long table instead of Deimos workbooks (CSV/XLSX)")
    p.add_argument("--sheet", help="sheet of --long-table")
    p.add_argument("--temperature-regex", help="regex giving the temperature id from the path, e.g. '(T\\d+)'")
    p.add_argument("--celsius", help="temperature → °C mapping (CSV: temperature,celsius)")
    p.add_argument("-o", "--out", help="output folder")
    p.add_argument("--fc", type=float, help="|log2 FC| threshold (default 0.6)")
    p.add_argument("--logp", type=float, help="-log10(p) threshold (default 2)")
    p.add_argument("--rule", choices=["strict", "r_legacy"], help="call rule (default strict)")
    p.add_argument("--min-temperatures", type=float, help="integer, or fraction of the series if < 1")
    p.add_argument("--min-streak", type=float, help="integer, or fraction of the series if < 1")
    p.add_argument("--no-plots", action="store_true", help="tables only")
    p.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                   help="override any key, e.g. --set call.min_consecutive=3 (repeatable)")
    p.add_argument("-v", "--verbose", action="store_true")
    return ap


def _apply_cli(cfg: dict, a: argparse.Namespace) -> dict:
    def put(attr, key):
        v = getattr(a, attr, None)
        if v not in (None, False):
            set_dotted(cfg, key, v)

    put("celsius", "temperatures_celsius")
    put("temperature_regex", "input.temperature_regex")
    put("out", "output.out_dir")
    put("fc", "call.fc_cutoff")
    put("logp", "call.logp_cutoff")
    put("rule", "call.rule")
    put("min_temperatures", "filter.min_temperatures")
    put("min_streak", "filter.min_streak")
    if a.deimos_root:
        cfg["input"].update(mode="deimos", deimos_root=a.deimos_root, manifest=None)
    if a.manifest:
        cfg["input"].update(mode="deimos", manifest=a.manifest)
    if a.long_table:
        cfg["input"].update(mode="long-table", long_table=a.long_table)
        if a.sheet:
            cfg["input"]["long_table_sheet"] = a.sheet
    if a.no_plots:
        for k in ("volcano", "profiles", "plot3d", "gif"):
            set_dotted(cfg, f"plots.{k}", False)
    for item in a.set or []:
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
        out = pipeline.analyze(cfg)
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        log.error("%s", e)
        return 1
    res = out["result"]
    print("\nSummary")
    print(res.summary.to_string(index=False))
    print(f"\nResults in {Path(cfg['output']['out_dir']).resolve()}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
