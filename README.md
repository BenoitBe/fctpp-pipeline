<p align="center">
  <img src="olympus_icon.svg" width="90" alt="Olympus">
</p>

# Olympus — fold-change Thermal Proteome Profiling

[![tests](https://github.com/BenoitBe/olympus-tpp-pipeline/actions/workflows/tests.yml/badge.svg)](https://github.com/BenoitBe/olympus-tpp-pipeline/actions/workflows/tests.yml)

Python port of the R package **FC.TPP**, named after Olympus Mons, the largest volcano
of the solar system: a nod to the temperature series, and its flank has the shape of
a melting curve (see the icon). Thermal-shift candidates are called from
per-temperature fold changes and p-values: a protein is a hit when it passes
|log2FC| and p-value thresholds at **several adjacent temperatures** with the same
direction.

Olympus starts from the **X result workbooks produced by
[Deimos](https://github.com/BenoitBe/deimos-proteomics-pipeline)**, one Deimos run
per temperature point (X = number of temperatures, any value). Deimos does the
per-temperature statistics (filtering, imputation, limma/eBayes); Olympus only
contains the TPP layer, so the platform keeps a single statistical engine.

```
Deimos × X (one run per temperature, run as usual)
   T1/ProteomicAnalysis_Results.xlsx … TX/ProteomicAnalysis_Results.xlsx
                         │  olympus analyze
                         ▼
   Olympus_Results.xlsx · volcano_plots.pdf · TPP_plots.pdf · <contrast>_plot3D.html
```

---

## Installation

```bash
git clone https://github.com/BenoitBe/olympus-tpp-pipeline.git
cd olympus-tpp-pipeline
pip install -e .            # or: pip install -r requirements.txt && pip install -e .
pip install adjustText      # optional: non-overlapping labels on volcano plots
```

Python ≥ 3.10. Deimos is not a dependency: only its result workbooks are read
(tested with Deimos commits `4d49d24`, 27/09/2026, and `e160463`, 17/11/2026).

---

## Quick start

```bash
olympus analyze --deimos-root deimos_results -o olympus_output
olympus analyze --deimos-root deimos_results --celsius temperatures.csv -o olympus_output
olympus init-config -o olympus_config.yaml && olympus analyze -c olympus_config.yaml   # YAML
```

Try it without data:

```bash
python example_data/make_synthetic.py            # 12 Deimos-like workbooks with planted shifts
olympus analyze --deimos-root example_data/deimos_results --celsius example_data/temperatures.csv -o example_output
```

---

## Inputs

**The X Deimos result workbooks.** One `ProteomicAnalysis_Results.xlsx` per
temperature, found under `--deimos-root`:

```
deimos_results/T1/ProteomicAnalysis_Results.xlsx                  # Deimos out_dir = one folder per temperature
deimos_results/T1/proteomics_output/ProteomicAnalysis_Results.xlsx
deimos_results/T1_ProteomicAnalysis_Results.xlsx                  # renamed copies in one folder
```

- **Temperature id.** It is taken from the first folder under the root (or from
  the file name). For other names, `--temperature-regex '(T\d+)'` takes it from the
  path, e.g. `TPP_T3_output/...` gives `T3`. A manifest also works:
  `--manifest runs.csv`, with columns `temperature,results_path[,temperature_c]`.
- **What is read.** From each workbook, the `Differential_Expression` sheet: `name`,
  `Genes`, `num_NAs`, and `<contrast>_diff`, `_p.val`, `_p.adj` for every contrast.
  The `Methods` sheet is read as well.
- **Order and temperature values.** Temperatures are ordered naturally (T1 < T2 < …
  < T10). You can also give their values in °C (`--celsius`, CSV
  `temperature,celsius`); they are then used for ordering and axes.

**Checklist for the X Deimos runs** (checked by Olympus, with a warning if not met):

| Point | Why | What Olympus does |
|---|---|---|
| Same condition order and same control at every temperature (e.g. `control_condition: Ctl`) | identical contrast names | a contrast reversed at one temperature is flipped back |
| Same processing: explicit `impute_method` (not `auto`, the Deimos default), same FDR mode | comparable temperatures | differences reported in `Methods` |
| Same primary statistic: a `report.pr_matrix.tsv` for **every** temperature or for none | Deimos uses `primary_statistic: dpc` (imputation-free) only when it finds a usable pr_matrix, and falls back to imputation otherwise | a mix of DPC and imputation is reported in `Methods` |
| One design per temperature, with that temperature's labels (`Sple-1B` at T1, `Sple-2B` at T2…) | exact label matching in Deimos | — |
| DIA-NN matrices with **6 annotation columns** (DIA-NN ≥ 2) | Deimos reads `columns[:6]` as annotations; an old export with 4 columns makes Deimos stop, reporting two samples as missing | — |

---

## Method

1. **Presence filter**: `nb_temp ≥ 6` temperatures where Deimos kept the protein, and
   `max_consecutive ≥ 6` adjacent temperatures (R `tidy_tpp`, same defaults).
   These defaults were set for 12 points. With another X, give numbers
   (`--min-temperatures 4 --min-streak 4`) or fractions of the series (`0.5` = half,
   rounded up). A filter larger than X stops the run with an explicit message.
2. **Per temperature and contrast**: `logpval = -log10(p)`. The p-value is `p.val` by
   default (as in R); `call.p_column: p.adj` switches to the Deimos BH value.
3. **Call** (`<contrast>_threshold = Over/Under`):
   - `strict` (default): at least `min_consecutive` (2) **adjacent** temperatures of
     the series, each with |log2FC| ≥ 0.6 and -log10 p ≥ 2, with the **same sign**.
   - `r_legacy`: exact R behaviour. Two *successive rows* are enough, even when a
     temperature is missing between them, and the sign is not checked.
4. **Flags on each hit** (sheet `Hits`):
   - `shift_class`:
     - `thermal`;
     - `abundance-like`: the difference is already there at the lowest temperature,
       so it is expression or solubility rather than stability;
     - `undetermined`: the protein is absent at the reference temperature;
     - `mixed-sign`: `r_legacy` only.
   - `imputed_in_run`: Deimos imputed at least one value of the protein (in any
     condition) at a temperature of the calling run.
   - `run_start`–`run_end`, `direction`, `n_sig_temps`, `max_abs_diff`, `max_logp`.
   - `over_strict` / `over_r_legacy`: the result of each rule, for comparison.

> **Multiple testing.** p-values are thresholded temperature by temperature. The
> consecutive-temperature rule limits false positives but does not calibrate them
> (adjacent temperatures come from the same lysates and are correlated). Treat hits
> as candidates, check them on the raw profiles, and read the flags.

---

## Outputs (`olympus_output/`)

| File | Content |
|---|---|
| `Olympus_Results.xlsx` | Sheets:<br>• `Methods`, `Summary`<br>• `Hits`, `TPP_data` (R layout), `Stats_proteins`, `Calls`<br>• `FC_<contrast>` / `logP_<contrast>`: protein × temperature matrices<br>• `Inputs`: path, SHA-256 and Deimos parameters for each temperature<br>• `Deimos_Methods`, `Temperatures`, `Parameters` |
| `TPP_data.csv`, `TPP_hits.csv` | Same tables, flat. |
| `volcano_plots.pdf` | One page per contrast, one panel per temperature (`volcano_layout: pages` gives the R layout). Hits are circled. |
| `TPP_plots.pdf` | FC and -log10 p profiles: all proteins in grey with hits highlighted, then one panel per hit. |
| `<contrast>_plot3D.html` | Interactive 3D view (temperature × log2FC × -log10 p) with threshold planes. `plotly_js: embed` (default) works offline. |
| `<contrast>_rotation.gif` | Optional (`plots.gif: true`), rendered with matplotlib. |
| `last_config.yaml` | Full configuration of the run. |

Consistency checks on the inputs:
- if a Deimos parameter (imputation, filtering, statistic…) differs between
  temperatures, you get a warning and a row in `Methods`;
- a contrast that comes out reversed at one temperature (`Ctl_vs_C1` instead of
  `C1_vs_Ctl`) is flipped back automatically.

---

## From the R package

| R (FC.TPP) | Python (Olympus) |
|---|---|
| `process_tsv`, `format_tsv`, `find_design`, `sumXP`, `process_and_export` (DEP: filter, MinProb, limma, fdrtool) | **Deimos**, one run per temperature (outside this package) |
| — | `inputs.discover_deimos_runs` / `load_deimos_runs`: reading the X workbooks |
| `tidy_tpp` | `tpp.run_tpp` (`rule: r_legacy` = identical calls) |
| `magma` | `plots.volcano_pdf` |
| `plot_and_3D` (ggplot2, plotly, webshot2, magick) | `plots.profiles_pdf`, `plot3d.plot3d_html`, `plot3d.rotation_gif` |
| `TPP_results.xlsx` (overwritten by each function) | `Olympus_Results.xlsx`, written once |

The R audit, with the list of changes and their justification, is in
[`docs/AUDIT_FC.TPP_R.md`](docs/AUDIT_FC.TPP_R.md) (in French).

---

## Validation

On the 12-temperature *S. aureus* dataset of the R project:

| Check | Result |
|---|---|
| Calls, `r_legacy` rule, on the R `TPP_data` sheet | identical: 23,792 rows, 3 contrasts, `_logpval` and `_threshold` |
| Presence statistics recomputed from the raw `T*.tsv` (cRAP, DEP `filter_missval`, label matching) | identical: 2,215 proteins, `nb_temp` and `max_consecutive`, 2,041 retained |
| Deimos vs DEP log2FC, rows without imputation | identical (max difference 2 × 10⁻¹⁴, 22,423 rows) |
| Deimos vs DEP log2FC, rows with imputation (QRILC vs MinProb) | median \|Δ\| = 0.29 (1,369 rows) |
| Hits | Both R hits found again (C1_vs_Ctl, T9–T10, up). One extra hit in C2_vs_Ctl, flagged `imputed_in_run` + `undetermined` (absent at T1). |

---

## Tests

```bash
pytest                                        # unit + CLI tests on synthetic data
OLYMPUS_R_PROJECT=/path/FC.TPP pytest -m realdata               # parity with the R outputs
OLYMPUS_DEIMOS_RESULTS=/path/deimos_results OLYMPUS_R_PROJECT=/path/FC.TPP pytest -m deimos   # real Deimos outputs
```

Real data are not versioned. The two marked suites are skipped when their variables
are not set.

---

## Limitations and ideas

- No FDR across temperatures. Options: permutation of the condition labels within
  replicates, or a curve model (NPARC-type, Childs et al. 2019).
- Per-temperature statistics cannot separate abundance from stability. A limma model
  with a temperature × condition interaction, fitted on all temperatures together,
  would. For now the `abundance-like` flag is a heuristic.
- Deimos side: reading `columns[:6]` as annotations and the `Sple-…_` regex are
  hard-coded (see the audit).

---

## Ecosystem

| Tool | Role |
|---|---|
| [Deimos](https://github.com/BenoitBe/deimos-proteomics-pipeline) | DIA quantitative pipeline (per-temperature statistics) |
| **Olympus** (this repo) | Thermal proteome profiling on top of Deimos (formerly FC.TPP, R) |
| [Phobos](https://github.com/BenoitBe/phobos-peptidomics-pipeline) | DDA/PEAKS peptidomics |
| [Pathfinder](https://github.com/BenoitBe/Pathfinder) | Cross-run synthesis |

## Authors and license

Benoit Bernay, Anthony Sleiman, Alexandre Mahé, with contributions from Aurélie
Budin-Verneuil (Université de Caen Normandie; funded by Région Normandie). Authors
and license are those of the R package: Artistic License 2.0, see [LICENSE](LICENSE).
See [CITATION.cff](CITATION.cff) to cite the software.
