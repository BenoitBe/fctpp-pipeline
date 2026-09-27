# FC.TPP — fold-change Thermal Proteome Profiling (Python)

Python port of the R package **FC.TPP**. Thermal-shift candidates are called from
per-temperature fold changes and p-values: a protein is a hit when it passes
|log2FC| and p-value thresholds at **several adjacent temperatures** with the same
direction.

The per-temperature statistics (filtering, imputation, limma/eBayes) are **no longer
computed in this package**. They come from
[Deimos](https://github.com/BenoitBe/deimos-proteomics-pipeline), run once per
temperature, so the platform keeps a single statistical engine. FC.TPP only
contains the TPP-specific layer.

```
T1.tsv … T12.tsv (DIA-NN)        fctpp prepare        deimos_runs/T1 … T12/
+ experimental design     ───────────────────────▶  pg_matrix (Deimos layout) + design + YAML
                                                              │  fctpp run-deimos (Deimos × 12)
                                                              ▼
                                                     ProteomicAnalysis_Results.xlsx × 12
                                                              │  fctpp analyze
                                                              ▼
                             FCTPP_Results.xlsx · volcano_plots.pdf · TPP_plots.pdf · *_plot3D.html
```

---

## Installation

```bash
git clone https://github.com/BenoitBe/fctpp-pipeline.git
cd fctpp-pipeline
pip install -e .            # or: pip install -r requirements.txt && pip install -e .
pip install adjustText      # optional: non-overlapping labels on volcano plots
```

Deimos is needed for the `run-deimos` step only (Python ≥ 3.10; tested with Deimos
commit `4d49d24`):

```bash
git clone https://github.com/BenoitBe/deimos-proteomics-pipeline.git
pip install -r deimos-proteomics-pipeline/requirements.txt
```

---

## Quick start

```bash
fctpp init-config -o fctpp_config.yaml       # annotated template, edit the paths

# 1. one Deimos input set per temperature
fctpp prepare --raw-dir TPP_data --design TPP_data/experimental_design.xlsx --runs deimos_runs

# 2. run Deimos on every temperature (≈ 30 s each; --jobs to parallelise)
fctpp run-deimos --runs deimos_runs --deimos-dir ../deimos-proteomics-pipeline --jobs 2

# 3. TPP calls, Excel report and figures
fctpp analyze --deimos-root deimos_runs -o fctpp_output

# or everything at once
fctpp all -c fctpp_config.yaml
```

If the Deimos runs already exist (run by hand, other settings…), step 3 is enough:
`fctpp analyze` looks for `<root>/<temperature>/**/ProteomicAnalysis_Results.xlsx`.
A manifest (`--manifest runs.csv`: `temperature,results_path[,temperature_c]`)
works too.

Try it without data:

```bash
python example_data/make_synthetic.py            # 12 Deimos-like workbooks with planted shifts
fctpp analyze --deimos-root example_data/deimos_runs --celsius example_data/temperatures.csv -o example_output
```

---

## Inputs

| File | Description |
|---|---|
| `TPP_data/T1.tsv … Tn.tsv` | One DIA-NN `pg_matrix` per temperature. Old exports (4 annotation columns) and DIA-NN ≥ 2 (6 columns) are both accepted. |
| Design (`.xlsx` / `.csv`) | `label;condition;replicate`, see below. |
| Temperatures in °C (optional) | `temperature,celsius` CSV or a dict in the YAML. Used for ordering and axes. |

**Design, two formats**

- *Template* (R FC.TPP format): the labels of **one** temperature (`Sple.1B`, `Sple.1C`…).
  The other temperatures are matched on the suffix that remains after removing the
  longest common prefix (`B`, `C`… `J`). This is what `DEP::make_se()` did silently in R.
  The match must be one-to-one, otherwise the run stops.
- *Long* (recommended): a `temperature` column plus explicit labels for every
  temperature.

Every match is written to `deimos_runs/sample_mapping.csv`.

`fctpp prepare` also:
- rewrites each matrix in the 6-annotation-column layout that Deimos expects (Deimos
  reads `columns[:6]` as annotations, so an old 4-column export would lose two runs);
- writes per-temperature designs with the exact labels Deimos derives from the run
  names (`(Sple-[^_]+)_`), in the same order at every temperature;
- writes one Deimos YAML per temperature, with no interactive module and a forced
  control (`control_condition: Ctl`), so that the contrasts have the same names and
  orientation at all temperatures.

---

## Method

1. **Presence filter**: `nb_temp ≥ 6` temperatures where Deimos kept the protein, and
   `max_consecutive ≥ 6` adjacent temperatures (R `tidy_tpp`, same defaults).
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
   - `imputed_in_run`: Deimos imputed at least one value at a temperature of the calling run.
   - `run_start`–`run_end`, `direction`, `n_sig_temps`, `max_abs_diff`, `max_logp`.
   - `over_strict` / `over_r_legacy`: the result of each rule, for comparison.

> **Multiple testing.** p-values are thresholded temperature by temperature. The
> consecutive-temperature rule limits false positives but does not calibrate them
> (adjacent temperatures come from the same lysates and are correlated). Treat hits
> as candidates, check them on the raw profiles, and read the flags.

---

## Outputs (`fctpp_output/`)

| File | Content |
|---|---|
| `FCTPP_Results.xlsx` | Sheets:<br>• `Methods`, `Summary`<br>• `Hits`, `TPP_data` (R layout), `Stats_proteins`, `Calls`<br>• `FC_<contrast>` / `logP_<contrast>`: protein × temperature matrices<br>• `Inputs`: path, SHA-256 and Deimos parameters for each temperature<br>• `Deimos_Methods`, `Temperatures`, `Parameters` |
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

| R (FC.TPP) | Python |
|---|---|
| `process_tsv`, `format_tsv`, `find_design`, `sumXP` | `fctpp prepare` (`prepare.py`, `design.py`) |
| `process_and_export` (DEP: filter, MinProb, limma, fdrtool) | **Deimos** (one run per temperature) |
| `tidy_tpp` | `tpp.run_tpp` (`rule: r_legacy` = identical calls) |
| `magma` | `plots.volcano_pdf` |
| `plot_and_3D` (ggplot2, plotly, webshot2, magick) | `plots.profiles_pdf`, `plot3d.plot3d_html`, `plot3d.rotation_gif` |
| `TPP_results.xlsx` (overwritten by each function) | `FCTPP_Results.xlsx`, written once |

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
FCTPP_R_PROJECT=/path/FC.TPP pytest -m realdata               # parity with the R outputs
FCTPP_DEIMOS_DIR=/path/deimos FCTPP_R_PROJECT=/path/FC.TPP pytest -m deimos   # real Deimos, slow
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
- Deimos' sample regex (`Sple-…_`) is hard-coded on the Deimos side.
  `prepare.sample_regex` has to stay in sync with it.

---

## Ecosystem

| Tool | Role |
|---|---|
| [Deimos](https://github.com/BenoitBe/deimos-proteomics-pipeline) | DIA quantitative pipeline (per-temperature statistics) |
| **FC.TPP** (this repo) | Thermal proteome profiling on top of Deimos |
| [Phobos](https://github.com/BenoitBe/phobos-peptidomics-pipeline) | DDA/PEAKS peptidomics |
| [Pathfinder](https://github.com/BenoitBe/Pathfinder) | Cross-run synthesis |

## Authors and license

Benoit Bernay, Anthony Sleiman, Alexandre Mahé, with contributions from Aurélie
Budin-Verneuil (Université de Caen Normandie; funded by Région Normandie). Authors
and license are those of the R package: Artistic License 2.0, see [LICENSE](LICENSE).
See [CITATION.cff](CITATION.cff) to cite the software.
