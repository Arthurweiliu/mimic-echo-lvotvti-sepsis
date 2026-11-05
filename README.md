# LVOT-VTI and LVEF in sepsis — analysis code

Analysis code for the manuscript:

> **LVOT-VTI and LVEF in sepsis: non-linear associations with 90-day mortality and the limited value of derived cardiac output**

The study uses MIMIC-IV v3.1 and MIMIC-IV-ECHO v1.0.1. **No patient-level data are included in this repository**; the MIMIC databases are available to credentialed users under a data use agreement (https://physionet.org/content/mimiciv/ and https://physionet.org/content/mimic-iv-echo/). This repository contains the code, the protocol, and the aggregate results needed to reproduce every number in the paper.

---

## 1. What the analysis does

1. **Cohort.** First ICU stay per patient with Sepsis-3, a structured LVOT-VTI measurement in the index transthoracic echocardiogram performed 0–72 h after ICU admission, alive at the index study, with left truncation at the study (n = 4,816; 1,573 deaths).
2. **Primary association.** VTI < 17 vs >= 17 cm with 90-day mortality, multiple imputation (m = 20) of covariates, Cox model with entry at the index study. Co-primary 72-h landmark; 7-day time split.
3. **Head-to-head.** VTI, LVEF and echocardiography-derived cardiac index (eCI) from the same study, modelled with 3-df natural cubic splines. Every competing model is fitted inside each imputation and predicted risks are averaged before discrimination/calibration are computed (H1 n = 4,336; H2 n = 3,815).
4. **Incremental value and validation.** delta AUC, NRI, IDI, decision curves with a paired bootstrap; bootstrap optimism-corrected C-index; entry-conditional calibration.

## 2. Requirements

| Component | Version used |
|---|---|
| R | 4.5.2 |
| mice | 3.19.0 |
| survival | 3.8.3 |
| data.table | 1.18.6.1 |
| splines, stats, graphics | base R |
| Node.js (audit and Markdown→DOCX only) | 24.x, package `docx` |

Total runtime on one core is roughly 2–3 hours, dominated by the two bootstrap jobs (`68_mi_delta.R`, `73_optimism.R`, ~13 min each) and the six measure-specific imputations in `79_s3_mi.R` (~25 min).

## 3. Inputs

Point `MIMIC_ECHO_DIR` at a directory holding the derived tables exported from your own MIMIC-IV / MIMIC-IV-ECHO instance. The build steps that produce them are in `scripts/00_inputs/`.

| File | Produced by | Content |
|---|---|---|
| `sepsis3_icustays.csv` | `00_inputs/build_sepsis3_icustays.sql` | Sepsis-3 ICU stays, admission times |
| `cohort_covariates.csv` | `00_inputs/build_cohort_covariates.sql` | first-day covariates, outcomes |
| `linked_cohort.csv` | `00_inputs/build_cohort_link.mjs` | ICU stay ↔ echo study link with signed hours from admission |
| `h1_data.csv` | `00_inputs/build_h1.mjs` | one row per stay: covariates, LVEF, LVOT diameter, VTI, follow-up |
| `analysis_data.csv` | `00_inputs/merge_analysis.mjs` | echo measurements, body surface area |
| `echo_extra.csv` | `00_inputs/extract_echo_extra.mjs` | long-format echo measurements (resting heart rate) |
| `lac.csv` | user export of `labevents` itemid 50813 | serum lactate (chemistry) |

**Linkage note.** `build_cohort_link.mjs` attaches to each ICU stay the echo study closest to ICU admission among studies within `[intime − 72 h, outtime + 72 h]`. The analysis then requires the *signed* interval 0–72 h, so a stay whose closest study preceded admission is excluded even if another study fell inside the window (79 stays, 1.6% of those with an in-window study). This is deliberate: it keeps a single, unambiguous time zero.

## 4. Layout

```
config.R                     paths (H1_WORK, MIMIC_ECHO_DIR, results/)
run_all.R                    runs every stage in order
scripts/00_inputs/           SQL / Node extracts that build the derived tables
scripts/01_cohort/           lactate supplement, flow, counts
scripts/02_primary/          primary MI model, time split, adjustment ladder, H1 subset
scripts/03_headtohead/       spline models, incremental value, bootstrap, optimism
scripts/04_sensitivity/      residual sensitivity, functional form, eCO alternative
scripts/05_tables_figures/   Table 1, Table 3-4, Figures 1-5
scripts/06_audit/            number audit, table check, Markdown → DOCX
results/                     aggregate outputs (no patient-level rows)
figures/                     final 600-dpi LZW TIFF figures
docs/                        protocol, cover letter, response to editorial points
```

## 5. Running it

```bash
export H1_WORK=/path/to/this/repository
export MIMIC_ECHO_DIR=/path/to/derived/tables
cd "$H1_WORK"

Rscript scripts/01_cohort/51_lactate_rebuild_v2.R     # -> results/lactate_v2.csv
Rscript scripts/02_primary/52_mi_primary_v2.R         # -> results/TableS_missing_data.csv, mi_comp_*.rds
...
```

or simply `Rscript run_all.R`. Each script starts with `source("config.R")` and therefore must be launched **from the repository root**; outputs go to `results/`.

Two long jobs take an environment variable for the number of bootstrap replicates:

```bash
SET=H1 BOOT_REPS=200 Rscript scripts/03_headtohead/68_mi_delta.R
SET=H2 BOOT_REPS=200 Rscript scripts/03_headtohead/68_mi_delta.R
SET=H1 BOOT_REPS=200 Rscript scripts/03_headtohead/73_optimism.R
SET=H2 BOOT_REPS=200 Rscript scripts/03_headtohead/73_optimism.R
```

All scripts set `seed = 20261105`, so repeated runs reproduce the same imputations and bootstrap draws.

## 6. What to check after a run

| Quantity | Value in the paper | Where |
|---|---|---|
| Primary cohort | 4,816 patients, 1,573 deaths | `01_cohort/66_counts.R` output |
| VTI < 17 vs >= 17, MI | 1.32 (1.19–1.47), FMI 0.004 | `02_primary/52_mi_primary_v2.R` |
| Per +5 cm | 0.93 (0.90–0.97) | same |
| Days 0–7 / day-7 landmark | 1.50 (1.26–1.78) / 1.34 (1.17–1.53) | `02_primary/55_ph_mi.R` |
| H1 combined delta AUC | +0.0128 (0.0084–0.0197) | `03_headtohead/67`, `68` |
| H2 all-spline delta AUC | +0.0150 (0.0097–0.0231) | same |
| Optimism-corrected C | 0.7015 (H1) / 0.7079 (H2) | `03_headtohead/73` |
| Per-SD VTI / LVEF / eCI | 0.779 / 1.004 / 1.204 | `03_headtohead/67` |
| Functional form (six measures) | Table S3 | `04_sensitivity/79_s3_mi.R` |
| Residual sensitivity | Table S4 | `04_sensitivity/78_s4_mi.R` |
| 28-day contrast | 1.31 (1.16-1.48) | `02_primary/82_compare_addons.R` |
| Alternative VTI cut-points | Table S6; 1.28-1.33 | `04_sensitivity/83_cutpoints_knots_d3.R` |
| 4-knot spline | 2-delta log-likelihood 48.0 (1 df) | same |
| Phenotype by VTI category | Table S7 | `05_tables_figures/84_phenotype.R` |
| Optimism-corrected calibration slope | 0.979 (H1) / 0.989 (H2) | `03_headtohead/86_calib_opt.R` |
| eCO alternative | Table S5 | `04_sensitivity/80_eco_mi.R` |

`06_audit/64_number_audit.mjs` re-checks that every hazard ratio quoted in the manuscript appears in one of these outputs.

## 7. Implementation notes for replicators

These two details cost us time and are easy to get wrong:

- **Concordance orientation.** With left-truncated `Surv(start, stop, event)` data, `survival::concordance(y ~ risk)` returns *1 − C* for a risk score in the survival 3.8.x build used here. The scripts use `reverse = TRUE`, which was verified against a manual Harrell C computed with the same left-truncation rule (0.6991 both).
- **Standardised effects for a log-transformed exposure.** A per-SD hazard ratio must use the standard deviation of the *modelled* term, i.e. `exp(beta × sd(log x))`, not `sd(x)`.

Additional conventions:

- Multiple imputation fills covariates only; **VTI, LVEF and eCI are never imputed**. All models in a comparison use identical patients.
- Predicted 90-day risk is entry-conditional: `p = 1 − exp(−[H0(90) − H0(start)] × exp(lp))`, computed per imputation and averaged before AUC, C-index, Brier, calibration and decision curves.
- eCI is derived (`SV × heart rate / BSA`) and is labelled as such throughout; values whose inputs fell outside physiological ranges (VTI 5–45 cm, LVOT 1.0–3.2 cm, heart rate 25–220 bpm, BSA 0.9–3.0 m², SV 10–250 mL) were set to missing rather than interpreted.
- First-day lactate is a blood-gas value; when missing it is supplemented from serum chemistry (maximum within 24 h). Agreement between the two sources in patients with both is r = 0.91.

## 8. Licence and citation

Code released under the MIT licence (see `LICENSE`). If you use it, please cite the manuscript above. The MIMIC databases remain subject to their own data use agreements.
