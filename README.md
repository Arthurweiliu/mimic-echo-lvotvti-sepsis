# Analysis code — LVOT-VTI and 90-day mortality in sepsis (MIMIC-IV-ECHO)

Code accompanying the manuscript: *Association of Left Ventricular Outflow Tract Velocity–Time Integral with 90-Day Mortality in Sepsis: A Retrospective Cohort Study Using MIMIC-IV-ECHO*.

## Data (not included)
This study used the **credentialed** datasets **MIMIC-IV** and **MIMIC-IV-ECHO**, distributed by PhysioNet (https://physionet.org). Under the PhysioNet data use agreement, **no patient-level or derived patient-level data are included in this repository**. To reproduce the analysis, obtain credentialed access and download `structured-measurement.csv.gz` plus the MIMIC-IV tables, and build the `mimiciv_derived` concepts (https://github.com/MIT-LCP/mimic-code).

## Working-directory convention
Every script reads and writes **in one data directory**.

- Node scripts: `node <script> [data_directory]` — defaults to the script's own directory.
- R scripts: `Rscript <script>` — no positional argument; each script sets its working directory to **its own directory** (`--file=` convention), so place `scripts/` next to your data files (or run with the data directory as the process working directory).
- SQL scripts: `psql -d <db> -v out='<data_directory>/<file>.csv' -f <script>.sql`.

So the simplest workflow is to copy `scripts/` next to your data files (or pass the directory explicitly).

## Requirements
- R ≥ 4.5.2 (`survival`, `splines`, `mice`). If your packages live outside the default library, set `R_LIBS_USER` — or `H1_LIB` — to that path.
- Node.js ≥ 20 (CSV cohort-construction scripts; no npm dependencies).
- PostgreSQL with MIMIC-IV and the `mimiciv_derived` concepts (`sepsis3`, `icustay_detail`, `sofa`, `first_day_sofa`, `sapsii`, `oasis`, `apsiii`, `charlson`, `first_day_lab`, `first_day_vitalsign`, `first_day_bg`, `kdigo_stages`, `rrt`, `ventilation`, `norepinephrine_equivalent_dose`, `cardiac_marker`).

## Order of execution
1. `build_sepsis3_icustays.sql`, `build_cohort_covariates.sql`, `build_aki_outcomes.sql`, `get_anchor.sql` — extract the Sepsis-3 ICU census, clinical covariates/outcomes and anchor year groups from PostgreSQL.
1b. `build_af_flag.sql` — atrial arrhythmia flags per stay from the ICU Heart Rhythm record (`chartevents` item 220048) and from admission ICD codes (I48/427.3) → `af_icu.csv`, `af_icd.csv`.
2. `build_tte_pheno.mjs` — build the TTE phenotype wide table from `structured-measurement.csv.gz`.
3. `build_cohort_link.mjs` — link TTE studies to the Sepsis-3 ICU stays (±72 h) → `linked_cohort.csv`.
4. `merge_analysis.mjs` — merge the echo and covariate tables → `analysis_data.csv`.
5. `build_h1.mjs` — build the analysis dataset → `h1_data_allstays.csv` (one row per eligible ICU stay). The 90-day outcome follows protocol section 8: survival time is time to death, or administrative censoring at day 90 when no death is recorded within 90 days.
6. `build_congestion.mjs` — right-heart/venous congestion burden score → `congestion_analysis.csv` (input for Table S3).
7. `15_make_first_icu.R` — retain one record per patient (first ICU stay) → writes **`h1_data.csv`**, the analysis file read by every script below (plus an identical `h1_data_first.csv` alias).
8. `07_h1.R`, `08_h1_outputs.R`, `09_vti_cutoff.R`, `10_h1_combo.R`, `11_mi.R`, `12_audit.R`, `13_ph_sens.R`, `14_unit_of_analysis.R`, `16_flowchart.R`, `17_reviewer_checks.R`, `18_km_tdhr.R`, `19_selection.R`, `20_review2.R`, `21_primary_table_spline.R`, `22_review4.R`, `23_s8.R`, `24_lvef_vs_vti.R`, `25_forest.R`, `26_evalue.R`, `27_validity.R` — descriptive statistics, primary/secondary models, dose–response, cut-point sensitivity, combined flow–congestion phenotype, multiple imputation, proportional-hazards sensitivity, unit-of-analysis, selection analysis, E-value, measurement-validity analyses, and figures.
9. `extract_echo_extra.mjs` — stream `structured-measurement.csv.gz` and extract the measurements absent from the phenotype table: **echo-time heart rate** (`resting_hr`) and mitral/tricuspid pulmonary-pulmonary regurgitation severity → `echo_extra.csv` (required by scripts 34–38). Accepts an optional data-directory argument.
10. `28_item3.R`, `29_item3_90d.R`, `38_fdr_panel.R` — **landmark analysis (co-primary)**, decision-curve / NRI / IDI incremental value (Set D), and the multiplicity-controlled (Benjamini–Hochberg) effect-modification panel.
11. `32_modifier_panel.R`, `33_modifier_joint.R` — effect-modification panel (atrial fibrillation, aortic/mitral/tricuspid valve disease, tachycardia, image quality) and joint AF × tachycardia models.
12. `34_co_ci.R`, `35_hr_check.R` — derived **stroke volume, cardiac output and cardiac index**, and the **per-beat vs per-minute** decomposition (echo-time heart rate from `echo_extra.csv`; body surface area from `analysis_data.csv`).
13. `36_figure_thesis.R` — **Figure 3** (per-beat flow vs cardiac output; also writes `Item4_thesis_forest.csv`).
14. `37_audit_methods.R` — methodological audit: echo-HR provenance, LVOT-diameter plausibility, subset-selection diagnostics, the CO/HR cancellation test, and C-index consistency across analysis sets.

> **Derived-input note.** Scripts 34–38 read `echo_extra.csv` (produced by step 9). Scripts 28–38 classify patients into four analysis sets: **A** primary (n = 5,410), **B** landmark (n = 5,099), **C** haemodynamic (n = 4,810), **D** incremental (n = 5,416). Hazard ratios and C-indices are set-specific and must not be interchanged between sets (see `Frozen_numbers_v2.md`).

## Table and figure provenance
| Output | Script |
|---|---|
| Table 1 (baseline by VTI category, n = 5410) | `21_primary_table_spline.R` |
| Table 2 (primary, co-primary, dose–response) | `07_h1.R`, `21_primary_table_spline.R` |
| Table 3 (effect modification) | `27_validity.R`, `38_fdr_panel.R` |
| Table 4 (incremental / decision-curve value) | `28_item3.R`, `29_item3_90d.R` |
| Table S1 (cut-point scan) | `09_vti_cutoff.R`, `22_review4.R`, `41_review2_details.R` |
| Table S2 (multiple imputation) | `11_mi.R` |
| Table S3 (PH diagnostics, time-transform) | `12_audit.R`, `13_ph_sens.R` |
| Table S4 (measurement-validity analyses) | `27_validity.R` |
| Table S5 (quantitative bias analysis) | `26_evalue.R`, `27_validity.R` |
| Table S6 (derived haemodynamic variables) | `34_co_ci.R` |
| Table S7 (haemodynamic subset vs excluded, SMD) | `42_baseline_smd.R` |
| Table S8 (effect modification and multiplicity) | `27_validity.R`, `38_fdr_panel.R` |
| Table S9 (flow–congestion phenotype) | `10_h1_combo.R` |
| Table S10 (unit of analysis) | `14_unit_of_analysis.R` |
| Table S11 (selection, ascertainment, derivation sensitivity; incl. additional sensitivity analyses) | `40_reviewer_sensitivity.R`, `41_review2_details.R` |
| Table S12 (TTE vs no TTE, SMD) | `42_baseline_smd.R` |
| Table S13 (LVOT-VTI availability, SMD) | `42_baseline_smd.R` |
| Figure 1 (flow chart) | `16_flowchart.R` |
| Figure 2 (Kaplan–Meier) | `18_km_tdhr.R` |
| Figure 3 (per-beat vs per-minute flow) | `36_figure_thesis.R` |
| Figure 4 (decision curve, 28-day) | `28_item3.R` |
| Figure S1 (spline) | `21_primary_table_spline.R` |
| Figure S2 (decision curve, 90-day) | `29_item3_90d.R` |
| Figure S3 (forest) | `25_forest.R` |
| Figure S4 (cut-point scan) | `09_vti_cutoff.R` |
| Figure S5 (time-dependent HR) | `18_km_tdhr.R` |

The LVEF head-to-head analysis (`24_lvef_vs_vti.R`) remains in the repository but is no longer reported as a table.

Figures are written as 600 dpi PNG (`png(..., res=600)`) and converted to LZW TIFF at 600 dpi for submission. All annotation text is black; colour is used only for series lines.

Originally examined as an exploratory congestion analysis: `05_congestion.R`, `06_scan.R`.

## Licence
Released under the **MIT License** — see `LICENSE`. If a different licence is preferred, replace that file before publication.
