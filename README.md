# Analysis code — LVOT-VTI and 90-day mortality in sepsis (MIMIC-IV-ECHO)

Code accompanying the manuscript: *Association of Left Ventricular Outflow Tract Velocity–Time Integral with 90-Day Mortality in Sepsis: A Retrospective Cohort Study Using MIMIC-IV-ECHO*.

## Data (not included)
This study used the **credentialed** datasets **MIMIC-IV** and **MIMIC-IV-ECHO**, distributed by PhysioNet (https://physionet.org). Under the PhysioNet data use agreement, **no patient-level or derived patient-level data are included in this repository**. To reproduce the analysis, obtain credentialed access and download `structured-measurement.csv.gz` plus the MIMIC-IV tables.

## Requirements
- R ≥ 4.5 (packages: survival, splines, mice; a user library may be needed: set `R_LIBS_USER`)
- Node.js ≥ 20 (for the CSV cohort-construction scripts)
- A local PostgreSQL instance loaded with MIMIC-IV (derived tables: `sepsis3`, `icustay_detail`, `sofa`, `first_day_sofa`, `sapsii`, `oasis`, `apsiii`, `charlson`, `first_day_lab`, `first_day_vitalsign`, `first_day_bg`, `kdigo_stages`, `rrt`, `ventilation`, `norepinephrine_equivalent_dose`, `cardiac_marker`)

## Order of execution
1. `scripts/build_cohort_covariates.sql`, `build_aki_outcomes.sql`, `get_anchor.sql` — extract clinical covariates/outcomes from PostgreSQL (psql `\copy`).
2. `scripts/build_tte_pheno.mjs` — build the TTE phenotype wide table from `structured-measurement.csv.gz`.
3. `scripts/build_cohort_link.mjs` — link TTE studies to Sepsis-3 ICU stays (±72 h).
4. `scripts/build_h1.mjs` — build the analysis dataset (`h1_data.csv`).
5. `scripts/15_make_first_icu.R` — retain one record per patient (first ICU stay).
6. `scripts/07_h1.R`, `08_h1_outputs.R`, `09_vti_cutoff.R`, `10_h1_combo.R`, `11_mi.R`, `13_ph_sens.R`, `14_unit_of_analysis.R`, `16_flowchart.R`, `17_reviewer_checks.R`, `18_km_tdhr.R`, `19_selection.R`, `21_primary_table_spline.R`, `22_review4.R`, `23_s8.R`, `24_lvef_vs_vti.R` — descriptive statistics, primary/secondary models, dose–response, cut-point sensitivity, combined flow–congestion phenotype, multiple imputation, proportional-hazards sensitivity, temporal validation, selection analysis, and figures.

Originally examined as an exploratory congestion analysis: `05_congestion.R`, `06_scan.R`.

## Licence
Code is released for research use. Choose an explicit licence (e.g. MIT) before publication; a licence file is intentionally not committed here for the authors to select.
