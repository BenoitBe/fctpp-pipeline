# Changelog — FC.TPP (Python)

## [0.1.0] — 2026-09-27

First Python release. Port of the R package FC.TPP (0.0.0.9000, 2025), with the
DEP part delegated to Deimos.

### Added
- `fctpp prepare`: one Deimos input set per temperature from the raw DIA-NN
  matrices. Includes the 6-column annotation layout, the per-temperature designs,
  an explicit sample mapping (`sample_mapping.csv`) and a Deimos YAML per temperature.
- `fctpp run-deimos`: batch runs of Deimos. stdin is closed, there is one log per
  temperature, and `--jobs` runs temperatures in parallel.
- `fctpp analyze`: TPP layer on Deimos outputs, or on any long table (e.g. the R
  `TPP_data` sheet).
- `strict` call rule (adjacent temperatures, same sign). The `r_legacy` rule
  reproduces the R calls exactly.
- Hit flags: `shift_class` (thermal / abundance-like / undetermined / mixed-sign),
  `imputed_in_run`, the calling run (`run_start`–`run_end`), direction.
- Checks that every temperature was processed the same way in Deimos. A contrast
  that comes out reversed at one temperature is flipped back automatically.
- Temperatures in °C (optional) for ordering and axes. The series is no longer
  fixed to T1..T12.
- Single Excel workbook with Methods, Summary, Hits, TPP_data, Stats_proteins,
  Calls, FC/logP matrices and input provenance (SHA-256).
- Volcano grid per contrast, FC/-log10 p profiles, lighter interactive 3D HTML,
  and a GIF without webshot2/ImageMagick.
- Tests: unit and CLI tests on synthetic data, parity with the R outputs
  (`-m realdata`), end-to-end with Deimos (`-m deimos`). GitHub Actions CI.

### Changed (vs R)
- Filtering, imputation and limma are done by Deimos: QRILC by default instead of
  MinProb, and BH `p.adj` instead of the DEP/fdrtool q-value.
- Volcano thresholds follow the call thresholds by default (R used 0.3 against 0.6
  for the calls).
- Figure labels use gene names; the accession is used when there is no gene name.

### Fixed (vs R)
- The R Excel export overwrote itself from one function to the next (only
  `TPP_data` and `Stats_proteins` survived).
- The default design folder `tpp_data` did not match `TPP_data` on case-sensitive
  file systems.
- The misleading "proteins were imputed" message (it counted the proteins removed
  by the filter) is gone; filtering is documented in the Deimos Methods sheet.
