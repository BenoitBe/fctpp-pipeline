import re

import pandas as pd
import pytest
import yaml

from fctpp.design import canonical, clean_sample_name, match_samples, read_design
from fctpp.prepare import DEIMOS_ANNOTATION_LAYOUT, prepare_deimos_inputs, split_pg_matrix
from fctpp.temperatures import TemperatureSeries, natural_sorted


def test_natural_sort_and_celsius_order():
    assert natural_sorted(["T10", "T2", "T1", "T12"]) == ["T1", "T2", "T10", "T12"]
    ts = TemperatureSeries.from_ids(["a", "b", "c"], {"a": 50, "b": 37, "c": 44})
    assert ts.ids == ("b", "c", "a") and ts.reference == "b"
    with pytest.raises(ValueError):
        TemperatureSeries.from_ids(["a", "b"], {"a": 37})


def test_canonical_and_clean_names():
    assert canonical("Sple-1B") == canonical("Sple.1B") == "Sple.1B"
    col = r"E:\DATAS\2025\CBSA\04-TPP-B13\250429-Sple-12D_Slot1-26_1_513.d"
    assert clean_sample_name(col) == "Sple-12D"          # same as Deimos load_data
    assert clean_sample_name("no_match_here") == "no_match_here"


def _design():
    return pd.DataFrame({"label": ["Sple.1D", "Sple.1C", "Sple.1B"],
                         "condition": ["C2", "C1", "Ctl"], "replicate": ["1", "1", "1"]})


def test_match_exact_then_suffix():
    cols_t1 = [r"x\250429-Sple-1B_S", r"x\250429-Sple-1C_S", r"x\250429-Sple-1D_S"]
    m = match_samples(cols_t1, _design())
    assert (m["match"] == "exact").all()
    assert m["sample_label"].tolist() == ["Sple-1D", "Sple-1C", "Sple-1B"]
    cols_t10 = [c.replace("-1", "-10") for c in cols_t1]
    m = match_samples(cols_t10, _design())
    assert m["match"].iloc[0].startswith("suffix")
    assert m["sample_label"].tolist() == ["Sple-10D", "Sple-10C", "Sple-10B"]
    assert m["condition"].tolist() == ["C2", "C1", "Ctl"]


def test_match_fails_loudly():
    cols = [r"x\250429-Sple-2B_S", r"x\250429-Sple-2C_S", r"x\250429-Sple-2X_S"]
    with pytest.raises(ValueError, match="Suffix matching failed"):
        match_samples(cols, _design())
    with pytest.raises(ValueError, match="exact"):
        match_samples(cols, _design(), strategy="exact")


def test_read_design_variants(tmp_path):
    p = tmp_path / "d.csv"
    p.write_text("Label;Condition;Replicate\nS1;Ctl;1.0\nS2;C1;1\n")
    d = read_design(p)
    assert list(d.columns[:3]) == ["label", "condition", "replicate"]
    assert d["replicate"].tolist() == ["1", "1"]
    p.write_text("label,condition\nS1,Ctl\n")
    with pytest.raises(ValueError, match="missing design column"):
        read_design(p)


def test_split_pg_matrix_detects_annotations():
    df = pd.DataFrame({"Protein.Group": ["A"], "Protein.Names": ["A_X"], "Genes": ["a"],
                       "First.Protein.Description": ["d"], "run1": [1.0], "run2": [2.0]})
    ann, samples = split_pg_matrix(df)
    assert samples == ["run1", "run2"] and len(ann) == 4


@pytest.mark.parametrize("four", [True, False])
def test_prepare_writes_deimos_compatible_inputs(raw_dir, tmp_path, four):
    raw, design = raw_dir
    if not four:   # regenerate with DIA-NN 2.x annotation columns
        from fctpp.synthetic import make_raw_pg_matrices
        raw, design = make_raw_pg_matrices(tmp_path / "raw2", four_annotation_columns=False)
    out = tmp_path / "runs"
    man = prepare_deimos_inputs(raw, design, out, deimos_params={"control_condition": "Ctl"},
                                celsius={"T1": 37, "T2": 40, "T3": 44, "T4": 47})
    assert man["temperature"].tolist() == ["T1", "T2", "T3", "T4"]
    assert man["temperature_c"].tolist() == [37, 40, 44, 47]
    for t in man["temperature"]:
        tsv = pd.read_csv(out / t / "report.pg_matrix.tsv", sep="\t")
        # Deimos takes the first 6 columns as annotations
        assert list(tsv.columns[:6]) == DEIMOS_ANNOTATION_LAYOUT
        assert tsv.shape[1] == 6 + 9
        des = pd.read_csv(out / t / "ExperimentalDesign.csv", sep=";", dtype=str)
        # labels = what Deimos derives from the run names (same regex)
        derived = [re.search(r"(Sple-[^_]+)_", c).group(1) for c in tsv.columns[6:]]
        assert des["label"].tolist() == derived
        assert des["condition"].tolist() == ["C2", "C1", "Ctl"] * 3   # template order kept
        cfg = yaml.safe_load((out / t / "deimos_config.yaml").read_text())
        assert cfg["control_condition"] == "Ctl" and cfg["make_dashboard"] is False
        assert cfg["out_dir"].endswith(f"{t}/deimos_output")
    mapping = pd.read_csv(out / "sample_mapping.csv")
    assert len(mapping) == 4 * 9


def test_prepare_with_long_design(raw_dir, tmp_path):
    raw, _ = raw_dir
    rows = []
    for t in ["T1", "T2", "T3", "T4"]:
        n = t[1:]
        for L, c, r in zip("BCDEFGHIJ", ["Ctl", "C1", "C2"] * 3, "111222333"):
            rows.append({"temperature": t, "label": f"Sple-{n}{L}", "condition": c, "replicate": r})
    d = tmp_path / "long.csv"
    pd.DataFrame(rows).to_csv(d, sep=";", index=False)
    man = prepare_deimos_inputs(raw, d, tmp_path / "runs")
    assert (man["match"] == "exact").all()
