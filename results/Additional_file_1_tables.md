# Additional file 1. Supplementary tables

Supplementary tables for: *LVOT-VTI and LVEF in sepsis: non-linear associations with 90-day mortality and the limited value of derived cardiac output*.

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
|---|---|---:|---:|---|---|---:|
| MI, m = 20 | Blood gas only | 4,816 | 1,573 | 1.33 (1.19-1.48) | 2.6e-07 | 0.004 |
| MI, m = 20 | Blood gas + chemistry (primary) | 4,816 | 1,573 | 1.32 (1.19-1.47) | 2.8e-07 | 0.004 |
| MI, m = 20 | Blood gas + chemistry + nearest | 4,816 | 1,573 | 1.32 (1.19-1.47) | 3.4e-07 | 0.002 |
| Complete case | Blood gas only | 3,573 | 1,262 | 1.26 (1.12-1.42) | - | - |
| Complete case | Blood gas + chemistry | 3,902 | 1,343 | 1.29 (1.15-1.44) | - | - |
| Complete case | Blood gas + chemistry + nearest | 4,346 | 1,449 | 1.30 (1.16-1.45) | - | - |

*All models adjusted for age, sex, SOFA, Charlson index, lactate and creatinine, with left truncation at the index TTE. FMI, fraction of missing information. The primary analysis is multiple imputation (m = 20) in the 4,816-patient cohort; complete-case rows are shown for comparison.*

---

## Table S3. Linear versus non-linear functional form for every same-study measure


| Measure | Contrast | n | Deaths | Linear HR (95% CI) | p linear | Non-linearity p | Overall p | delta AIC (spline - linear) |
|---|---|---:|---:|---|---:|---:|---:|---:|
| LVOT-VTI | per 5 | 4,851 | 1,608 | 0.93 (0.90-0.97) | <0.001 | 0.00093 | 3.1e-08 | +10.0 |
| Stroke volume | per 10 | 4,579 | 1,508 | 0.97 (0.95-0.99) | 0.003 | 0.00019 | 2.3e-07 | +13.1 |
| Cardiac output | per 1 | 4,534 | 1,500 | 1.01 (0.99-1.03) | 0.440 | 0.0009 | 0.0024 | +10.0 |
| Cardiac index | per 1 | 4,307 | 1,413 | 1.06 (1.01-1.11) | 0.018 | 0.00039 | 0.00069 | +11.7 |
| Heart rate | per 10 | 4,803 | 1,598 | 1.08 (1.05-1.11) | <0.001 | 0.99 | 1.9e-07 | -4.0 |
| LVEF | per 10 | 4,369 | 1,436 | 0.97 (0.94-1.00) | 0.049 | 6.1e-12 | 5.4e-12 | +47.6 |

*Covariates imputed (m = 20) within each measure-specific cohort; the markers themselves are not imputed. The linear coefficient is fitted on the untransformed exposure (the quantity under test) and pooled by Rubin's rules; the non-linearity test compares a 3-df natural spline with a linear term on the same scale (log for the flow measures, untransformed for LVEF and heart rate) and is reported as the median p across the 20 imputations. A positive delta AIC means the spline fits better.*


---

## Table S4. Residual sensitivity analyses for the primary contrast (same imputation model)

| Analysis | n | Deaths | HR for VTI < 17 vs >= 17 cm (95% CI) |
|---|---:|---:|---|
| Primary (index TTE 0-72 h, first stay) | 4,851 | 1,608 | 1.32 (1.19-1.47) |
| All in-window stays (first stay not required) | 5,344 | 1,797 | 1.30 (1.18-1.43) |
| Index TTE within 24 h | 3,048 | 1,088 | 1.27 (1.12-1.44) |
| Index TTE within 48 h | 4,202 | 1,434 | 1.30 (1.17-1.46) |
| Excluding poor or suboptimal image quality | 3,258 | 1,082 | 1.42 (1.25-1.61) |
| Excluding moderate-or-worse aortic stenosis | 4,574 | 1,486 | 1.31 (1.17-1.46) |
| Excluding moderate or large pericardial effusion | 4,820 | 1,593 | 1.33 (1.19-1.47) |

**Note.** Every row is a separate multiple imputation (m = 20) in the stated subset with the same covariate set and imputation model as the primary analysis; the imaging-window rows restrict the index TTE further, and the last three rows exclude defined abnormalities. The primary row is repeated for reference.

---

## Table S5. Echocardiography-derived CO (eCO, no body surface area) instead of eCI

| Model | n | Deaths | AUC | delta AUC |
|---|---:|---:|---:|---:|
| Clinical | 4,054 | 1,326 | 0.7344 | - |
| VTI + LVEF splines | 4,054 | 1,326 | 0.7480 | +0.0136 |
| VTI + LVEF + eCO splines | 4,054 | 1,326 | 0.7484 | +0.0140 |

**Note.** eCO = stroke volume x echo-time heart rate (L/min) and does not require body surface area, so the analysis set is larger than H2 (4,054 vs 3,845 patients). Covariates imputed (m = 20); predicted risks averaged over imputations. eCO added +0.0004 beyond VTI + LVEF.

---

## Table S6. Alternative boundaries for the primary VTI contrast

| Definition | n low | n high | HR (95% CI) |
|---|---:|---:|---|
| VTI <17 cm vs rest | 1,344 | 3,507 | 1.325 (1.192-1.473) |
| VTI <=17 cm vs rest | 1,614 | 3,237 | 1.277 (1.152-1.416) |
| VTI <16.5 cm vs rest | 1,343 | 3,508 | 1.326 (1.192-1.474) |
| VTI <18 cm vs rest | 1,614 | 3,237 | 1.277 (1.152-1.416) |
| VTI <16 cm vs rest | 1,076 | 3,775 | 1.333 (1.192-1.490) |

**Note.** All five definitions were fitted in the same 20 imputations and the same patients; only the boundary moves. VTI is recorded to whole centimetres, so the "< 18" and "<= 17" definitions contain the same patients.

---

## Table S7. Haemodynamic and geometric phenotype by VTI category


| Variable | Overall | < 17 cm | 17-20 cm | > 20 cm | p | Missing, n |
|---|---|---|---|---|---:|---:|
| Heart rate, bpm | 86.4 (75.0-98.6) | 92.6 (80.7-106.4) | 87.5 (76.6-99.2) | 82.3 (71.6-93.3) | <0.001 | 3 |
| Mean arterial pressure, mmHg | 75.2 (69.6-82.6) | 75.9 (70.4-82.8) | 75.3 (70.4-82.5) | 74.7 (68.6-82.5) | <0.001 | 3 |
| Stroke volume, mL | 65.2 (50.3-81.7) | 42.5 (34.6-50.3) | 59.7 (53.4-68.4) | 79.7 (69.1-94.2) | <0.001 | 272 |
| Cardiac output, L/min | 5.50 (4.25-6.94) | 3.92 (3.11-4.76) | 5.24 (4.33-6.26) | 6.52 (5.46-7.95) | <0.001 | 275 |
| Estimated SVR, dyn s cm-5 | 1111 (867-1451) | 1546 (1263-2010) | 1176 (962-1433) | 923 (752-1123) | <0.001 | 275 |
| LV end-diastolic dimension, cm | 4.50 (4.00-5.10) | 4.60 (4.00-5.30) | 4.50 (3.90-5.00) | 4.50 (4.00-5.00) | <0.001 | 217 |
| Norepinephrine, ug/kg/min | 0.250 (0.120-0.525) | 0.301 (0.150-0.651) | 0.207 (0.100-0.501) | 0.201 (0.100-0.500) | <0.001 | 2192 |
| Lactate, mmol/L | 2.2 (1.4-3.8) | 2.8 (1.8-5.1) | 2.1 (1.4-3.5) | 1.9 (1.3-3.4) | <0.001 | 1244 |
| Urine output, first 24 h, mL | 1420 (795-2310) | 1278 (650-2121) | 1460 (875-2410) | 1482 (836-2350) | <0.001 | 130 |

*Median (IQR); Kruskal-Wallis test across the three categories. Cardiac output = stroke volume x echo heart rate / 1000; estimated SVR = 80 x mean arterial pressure / cardiac output. LV end-diastolic dimension was available only where the reporting template recorded it. Urine output is the first-24 h total from the derived covariates table.*

---

## Table S8. Cohort-definition sensitivity analysis

| Analysis | n | Deaths | HR for VTI < 17 vs >= 17 cm (95% CI) |
|---|---:|---:|---|
| Primary (first Sepsis-3 ICU stay) | 4,851 | 1,608 | 1.32 (1.19-1.47) |
| Restricted to the first ICU stay of any kind | 4,305 | 1,420 | 1.32 (1.18-1.48) |

**Note.** The first row is the primary analysis. The second implements a stricter unit-of-analysis rule: the patient's first ICU stay of any kind, which excludes 546 patients (11.3%) who entered through a later ICU stay after an earlier non-septic stay (90-day mortality 34.4% in those excluded versus 33.0% in the patients retained). Both models used the same covariates and multiple imputation (m = 20); death times were the precise hospital death timestamps where available (see Outcomes).
