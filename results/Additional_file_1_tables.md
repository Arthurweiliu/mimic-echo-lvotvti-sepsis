# Additional file 1. Supplementary tables

Supplementary tables for: *Beyond linear coefficients: LVOT-VTI and LVEF carry complementary non-linear prognostic information in sepsis, whereas echocardiography-derived cardiac output adds little*.

---

## Table S1. Selection by echocardiography: Sepsis-3 first ICU stays with and without a transthoracic echocardiogram

| Group | n | Age, years (median) | Female, % | SOFA (median) | 90-day mortality, % | 72-h mortality, % |
|---|---:|---:|---:|---:|---:|---:|
| No TTE | 13,476 | 66 | 41.4 | 3 | 20.5 | 6.0 |
| TTE without VTI | 4,867 | 67 | 41.3 | 3 | 32.5 | 6.9 |
| TTE with VTI | 7,253 | 68 | 44.8 | 3 | 29.9 | 4.5 |
| p | - | < 0.001 | < 0.001 | < 0.001 | < 0.001 | < 0.001 |

**Note.** Patient-level comparison of all Sepsis-3 first ICU stays in the source census. TTE denotes any transthoracic echocardiogram; VTI denotes a structured LVOT-VTI measurement. Age, female sex and mortality differed across groups, confirming selection by echocardiography.

---

## Table S2. Missing-data sensitivity: VTI < 17 vs >= 17 cm for 90-day mortality

| Analysis | Lactate source | n | Deaths | HR (95% CI) | p | FMI |
|---|---|---:|---:|---|---:|---:|
| MI, m = 20 | Blood gas only | 4,816 | 1,573 | 1.33 (1.19-1.48) | 2.6e-07 | 0.004 |
| MI, m = 20 | Blood gas + chemistry (primary) | 4,816 | 1,573 | 1.32 (1.19-1.47) | 2.8e-07 | 0.004 |
| MI, m = 20 | Blood gas + chemistry + nearest | 4,816 | 1,573 | 1.32 (1.19-1.47) | 3.4e-07 | 0.002 |
| Complete case | Blood gas only | 3,573 | 1,262 | 1.26 (1.12-1.42) | - | - |
| Complete case | Blood gas + chemistry | 3,902 | 1,343 | 1.29 (1.15-1.44) | - | - |
| Complete case | Blood gas + chemistry + nearest | 4,346 | 1,449 | 1.30 (1.16-1.45) | - | - |

**Note.** All models adjusted for age, sex, SOFA, Charlson index, lactate and creatinine, with left truncation at the index TTE. FMI, fraction of missing information. The primary analysis is multiple imputation (m = 20) in the 4,816-patient cohort. "Chemistry" denotes the maximum serum chemistry lactate within 24 h of ICU admission; "nearest" denotes the closest chemistry value within +/-12 h of that window.

---

## Table S3. Linear versus non-linear functional form for every same-study measure


| Measure | Contrast | n | Deaths | Linear HR (95% CI) | p linear | Non-linearity p | Overall p | delta AIC (spline - linear) |
|---|---|---:|---:|---|---:|---:|---:|---:|
| LVOT-VTI | per 5 | 4,816 | 1,573 | 0.93 (0.90-0.97) | <0.001 | 0.00085 | 3.8e-08 | +10.1 |
| Stroke volume | per 10 | 4,547 | 1,476 | 0.97 (0.95-0.99) | 0.006 | 0.00019 | 5.8e-07 | +13.2 |
| Cardiac output | per 1 | 4,502 | 1,468 | 1.01 (0.99-1.04) | 0.302 | 0.00057 | 0.0018 | +10.9 |
| Cardiac index | per 1 | 4,276 | 1,382 | 1.06 (1.02-1.11) | 0.008 | 0.00027 | 0.00037 | +12.5 |
| Heart rate | per 10 | 4,768 | 1,563 | 1.08 (1.05-1.11) | <0.001 | 0.97 | 2.1e-07 | -3.9 |
| LVEF | per 10 | 4,336 | 1,403 | 0.97 (0.94-1.00) | 0.065 | 1.9e-12 | 2.1e-12 | +50.0 |


**Note.** Covariates imputed (m = 20) within each measure-specific cohort; the markers themselves are not imputed. The linear coefficient is fitted on the untransformed exposure (the quantity under test) and pooled by Rubin's rules; the non-linearity test compares a 3-df natural spline with a linear term on the same scale (log for the flow measures, untransformed for LVEF and heart rate) and is reported as the median p across the 20 imputations. A positive delta AIC means the spline fits better.

---

## Table S4. Residual sensitivity analyses for the primary contrast (same imputation model)

| Analysis | n | Deaths | HR for VTI < 17 vs >= 17 cm (95% CI) |
|---|---:|---:|---|
| Primary (index TTE 0-72 h, first stay) | 4,816 | 1,573 | 1.32 (1.19-1.47) |
| All in-window stays (first stay not required) | 5,306 | 1,759 | 1.29 (1.17-1.43) |
| Index TTE within 24 h | 3,020 | 1,060 | 1.26 (1.11-1.44) |
| Index TTE within 48 h | 4,167 | 1,399 | 1.30 (1.17-1.46) |
| Excluding poor or suboptimal image quality | 3,234 | 1,058 | 1.42 (1.24-1.62) |
| Excluding moderate-or-worse aortic stenosis | 4,543 | 1,455 | 1.31 (1.17-1.46) |
| Excluding moderate or large pericardial effusion | 4,785 | 1,558 | 1.32 (1.19-1.47) |

**Note.** Every row is a separate multiple imputation (m = 20) in the stated subset with the same covariate set and imputation model as the primary analysis; the imaging-window rows restrict the index TTE further, and the last three rows exclude defined abnormalities. The primary row is repeated for reference.

---

## Table S5. Echocardiography-derived CO (eCO, no body surface area) instead of eCI

| Model | n | Deaths | AUC | delta AUC |
|---|---:|---:|---:|---:|
| Clinical | 4,023 | 1,295 | 0.7330 | - |
| VTI + LVEF splines | 4,023 | 1,295 | 0.7466 | +0.0135 |
| VTI + LVEF + eCO splines | 4,023 | 1,295 | 0.7471 | +0.0141 |

**Note.** eCO = stroke volume x echo-time heart rate (L/min) and does not require body surface area, so the analysis set is larger than H2 (4,023 vs 3,815 patients). Covariates imputed (m = 20); predicted risks averaged over imputations. eCO added +0.0006 beyond VTI + LVEF.
