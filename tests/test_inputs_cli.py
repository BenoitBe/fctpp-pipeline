import shutil

import openpyxl
import pandas as pd
import pytest

from olympus.cli import main
from olympus.inputs import (discover_deimos_runs, harmonise_contrasts, load_deimos_runs,
                          load_long_table, read_manifest)
from olympus.synthetic import CELSIUS, make_deimos_like_runs


def test_discover_and_load(deimos_root):
    runs = discover_deimos_runs(deimos_root)
    assert list(runs) == [f"T{i}" for i in range(1, 13)]
    li = load_deimos_runs(runs)
    assert li.contrasts == ["C2_vs_C1", "C2_vs_Ctl", "C1_vs_Ctl"]
    assert li.results["imputed"].dropna().isin([True, False]).all()   # VRAI/FAUX parsed
    assert len(li.provenance) == 12 and li.provenance["sha256"].str.len().eq(64).all()
    assert not li.warnings


def test_discovery_layouts(tmp_path):
    # flat folder of renamed copies, X = 5 temperatures, plus an Excel lock file
    flat = make_deimos_like_runs(tmp_path / "flat", n_proteins=40,
                                 temps=["T1", "T2", "T3", "T10", "T11"], layout="flat")
    (flat / "~$T1_ProteomicAnalysis_Results.xlsx").write_bytes(b"lock")
    assert list(discover_deimos_runs(flat)) == ["T1", "T2", "T3", "T10", "T11"]
    # folders named freely: the temperature id comes from a regex
    free = tmp_path / "free"
    for t in ("T1", "T2", "T12"):
        src = make_deimos_like_runs(tmp_path / f"src_{t}", n_proteins=20, temps=[t])
        dst = free / f"TPP_{t}_deimos" / "proteomics_output"
        dst.mkdir(parents=True)
        shutil.copy(src / t / "proteomics_output" / "ProteomicAnalysis_Results.xlsx", dst)
    runs = discover_deimos_runs(free, temperature_regex=r"_(T\d+)_")
    assert list(runs) == ["T1", "T2", "T12"]
    assert list(discover_deimos_runs(free)) == ["TPP_T1_deimos", "TPP_T2_deimos", "TPP_T12_deimos"]
    # two workbooks for one temperature (e.g. Deimos species split) -> explicit error
    extra = free / "TPP_T1_deimos" / "proteomics_output" / "STAA8"
    extra.mkdir()
    shutil.copy(runs["T1"], extra)
    with pytest.raises(ValueError, match="Two Deimos workbooks"):
        discover_deimos_runs(free, temperature_regex=r"_(T\d+)_")


def test_harmonise_flips_reversed_contrast():
    a = pd.DataFrame({"name": ["P"], "A_vs_B_diff": [1.0], "A_vs_B_p.val": [0.01]})
    b = pd.DataFrame({"name": ["P"], "B_vs_A_diff": [1.0], "B_vs_A_p.val": [0.01]})
    out, contrasts, warns = harmonise_contrasts({"T1": a, "T2": b})
    assert contrasts == ["A_vs_B"]
    assert out["T2"]["A_vs_B_diff"].iloc[0] == -1.0
    assert any("flipped" in w for w in warns)


def test_inconsistent_deimos_methods_are_reported(deimos_root, tmp_path):
    root = tmp_path / "runs"
    shutil.copytree(deimos_root, root)
    wb_path = root / "T5" / "proteomics_output" / "ProteomicAnalysis_Results.xlsx"
    wb = openpyxl.load_workbook(wb_path)
    ws = wb["Methods"]
    for row in ws.iter_rows():
        if row[0].value == "Imputation":
            row[1].value = "Mixed (MNAR+MAR)"
    wb.save(wb_path)
    li = load_deimos_runs(discover_deimos_runs(root))
    assert any("Imputation" in w and "differs" in w for w in li.warnings)


def test_manifest_and_long_table(deimos_root, tmp_path):
    runs = discover_deimos_runs(deimos_root)
    man = tmp_path / "runs.csv"
    pd.DataFrame({"temperature": list(runs), "results_path": [str(p) for p in runs.values()],
                  "temperature_c": [CELSIUS[t] for t in runs]}).to_csv(man, index=False)
    r, cel = read_manifest(man)
    assert list(r) == list(runs) and cel["T12"] == 70.0
    li = load_deimos_runs(r)
    long = tmp_path / "long.csv"
    li.results.to_csv(long, index=False)
    lt = load_long_table(long)
    assert lt.contrasts == li.contrasts and len(lt.results) == len(li.results)


def test_cli_analyze_end_to_end(deimos_root, tmp_path):
    out = tmp_path / "out"
    cel = tmp_path / "celsius.csv"
    pd.DataFrame({"temperature": list(CELSIUS), "celsius": list(CELSIUS.values())}).to_csv(cel, index=False)
    code = main(["analyze", "--deimos-root", str(deimos_root), "-o", str(out), "--celsius", str(cel),
                 "--set", "plots.gif=true", "--set", "plots.plotly_js=cdn"])
    assert code == 0
    for f in ["Olympus_Results.xlsx", "TPP_data.csv", "TPP_hits.csv", "volcano_plots.pdf",
              "TPP_plots.pdf", "C1_vs_Ctl_plot3D.html", "C1_vs_Ctl_rotation.gif", "last_config.yaml"]:
        assert (out / f).exists(), f
    wb = openpyxl.load_workbook(out / "Olympus_Results.xlsx", read_only=True)
    for s in ["Methods", "Summary", "Hits", "TPP_data", "Stats_proteins", "Calls",
              "FC_C1_vs_Ctl", "logP_C1_vs_Ctl", "Inputs", "Temperatures"]:
        assert s in wb.sheetnames, s
    data = pd.read_csv(out / "TPP_data.csv")
    assert "Temperature_C" in data.columns
    html = (out / "C1_vs_Ctl_plot3D.html").read_text(encoding="utf-8")
    assert len(html) < 2_000_000 and "gen0001" in html


def test_cli_rule_and_errors(deimos_root, tmp_path, capsys):
    assert main(["analyze", "--deimos-root", str(deimos_root), "-o", str(tmp_path / "o"),
                 "--rule", "r_legacy", "--no-plots"]) == 0
    hits = pd.read_csv(tmp_path / "o" / "TPP_hits.csv")
    assert {"P0010", "P0011"} <= set(hits.loc[hits["contrast"] == "C1_vs_Ctl", "name"])
    assert main(["analyze", "--deimos-root", str(tmp_path / "missing"), "-o", str(tmp_path / "o2")]) == 1


def test_cli_small_series_with_fractions(tmp_path):
    """X = 5 temperatures: the R defaults (6/6) are impossible and stop the run
    with an explicit message; fractions adapt the filter to the series."""
    root = make_deimos_like_runs(tmp_path / "x5", temps=[f"T{i}" for i in range(6, 11)])
    out = tmp_path / "o5"
    assert main(["analyze", "--deimos-root", str(root), "-o", str(out), "--no-plots"]) == 1
    assert main(["analyze", "--deimos-root", str(root), "-o", str(out), "--no-plots",
                 "--min-temperatures", "0.5", "--min-streak", "0.5"]) == 0
    hits = pd.read_csv(out / "TPP_hits.csv")
    assert "P0001" in set(hits["name"])


def test_init_config(tmp_path):
    p = tmp_path / "cfg.yaml"
    assert main(["init-config", "-o", str(p)]) == 0
    from olympus.config import load_config
    cfg = load_config(p)
    assert cfg["call"]["rule"] == "strict"
    assert cfg["input"]["workbook_glob"] == "**/*ProteomicAnalysis_Results*.xlsx"
    assert "prepare" not in cfg
