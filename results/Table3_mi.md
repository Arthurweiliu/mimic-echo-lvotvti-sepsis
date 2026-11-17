### Table 3. Head-to-head models in the multiply imputed cohort (natural cubic splines, 3 df)

| Model | C-index | AUC | delta AUC | Calibration slope | Intercept | Brier | NRI-cat | IDI |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **H1 (n = 4,369; 1,436 deaths)** | | | | | | | | |
| Clinical | 0.6948 | 0.7311 | - | 0.990 | -0.006 | 0.1889 | - | - |
| VTI spline | 0.7014 | 0.7389 | +0.0078 | 0.996 | -0.002 | 0.1872 | 0.026 | 0.0072 |
| LVEF spline | 0.7013 | 0.7386 | +0.0076 | 0.986 | -0.010 | 0.1866 | 0.028 | 0.0112 |
| VTI + LVEF splines | 0.7058 | 0.7440 | +0.0129 | 0.991 | -0.006 | 0.1856 | 0.037 | 0.0153 |
| **H2 (n = 3,845; 1,250 deaths)** | | | | | | | | |
| Clinical | 0.7008 | 0.7375 | - | 1.004 | 0.005 | 0.1867 | - | - |
| VTI spline | 0.7069 | 0.7450 | +0.0075 | 1.016 | 0.014 | 0.1850 | 0.021 | 0.0062 |
| LVEF spline | 0.7072 | 0.7451 | +0.0076 | 1.000 | 0.002 | 0.1844 | 0.026 | 0.0109 |
| VTI + LVEF splines | 0.7111 | 0.7497 | +0.0122 | 1.008 | 0.008 | 0.1835 | 0.038 | 0.0143 |
| eCI spline | 0.7021 | 0.7388 | +0.0013 | 1.004 | 0.004 | 0.1861 | 0.006 | 0.0023 |
| VTI + eCI splines | 0.7096 | 0.7484 | +0.0109 | 1.007 | 0.008 | 0.1839 | 0.032 | 0.0123 |
| LVEF + eCI splines | 0.7075 | 0.7453 | +0.0078 | 0.999 | 0.001 | 0.1842 | 0.026 | 0.0120 |
| All splines | 0.7130 | 0.7521 | +0.0146 | 1.001 | 0.004 | 0.1826 | 0.047 | 0.0189 |

*Predicted 90-day risks are the mean of the risks from the 20 imputations; discrimination and calibration are computed on that averaged risk. AUC, area under the receiver-operating-characteristic curve for 90-day mortality; C-index, Harrell concordance with left truncation; Brier score at 90 days; calibration slope and intercept from a logistic model of the outcome on the logit of the predicted risk; NRI-cat, three-level net reclassification improvement (<10%, 10-25%, >=25%) versus the clinical model; IDI, integrated discrimination improvement versus the clinical model. The clinical model contains age, sex, SOFA, Charlson index, lactate and creatinine.*
