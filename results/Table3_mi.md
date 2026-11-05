### Table 3. Head-to-head models in the multiply imputed cohort (natural cubic splines, 3 df)

| Model | C-index | AUC | delta AUC | Calibration slope | Intercept | Brier | NRI-cat | IDI |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **H1 (n = 4,336; 1,403 deaths)** | | | | | | | | |
| Clinical | 0.6934 | 0.7296 | - | 1.005 | 0.006 | 0.1884 | - | - |
| VTI spline | 0.6998 | 0.7372 | +0.0076 | 1.010 | 0.010 | 0.1867 | 0.029 | 0.0070 |
| LVEF spline | 0.7000 | 0.7372 | +0.0076 | 0.997 | 0.000 | 0.1861 | 0.036 | 0.0117 |
| VTI + LVEF splines | 0.7045 | 0.7425 | +0.0128 | 1.003 | 0.004 | 0.1851 | 0.043 | 0.0156 |
| **H2 (n = 3,815; 1,220 deaths)** | | | | | | | | |
| Clinical | 0.6991 | 0.7355 | - | 1.019 | 0.018 | 0.1863 | - | - |
| VTI spline | 0.7051 | 0.7430 | +0.0075 | 1.031 | 0.028 | 0.1847 | 0.019 | 0.0060 |
| LVEF spline | 0.7056 | 0.7431 | +0.0076 | 1.012 | 0.012 | 0.1840 | 0.025 | 0.0113 |
| VTI + LVEF splines | 0.7095 | 0.7478 | +0.0123 | 1.021 | 0.020 | 0.1831 | 0.037 | 0.0146 |
| eCI spline | 0.7004 | 0.7368 | +0.0014 | 1.018 | 0.017 | 0.1857 | -0.001 | 0.0024 |
| VTI + eCI splines | 0.7080 | 0.7466 | +0.0112 | 1.021 | 0.021 | 0.1834 | 0.027 | 0.0126 |
| LVEF + eCI splines | 0.7060 | 0.7435 | +0.0081 | 1.011 | 0.011 | 0.1837 | 0.020 | 0.0125 |
| All splines | 0.7115 | 0.7504 | +0.0150 | 1.014 | 0.015 | 0.1822 | 0.046 | 0.0196 |

*Predicted 90-day risks are the mean of the risks from the 20 imputations; discrimination and calibration are computed on that averaged risk. AUC, area under the receiver-operating-characteristic curve for 90-day mortality; C-index, Harrell concordance with left truncation; Brier score at 90 days; calibration slope and intercept from a logistic model of the outcome on the logit of the predicted risk; NRI-cat, three-level net reclassification improvement (<10%, 10-25%, >=25%) versus the clinical model; IDI, integrated discrimination improvement versus the clinical model. The clinical model contains age, sex, SOFA, Charlson index, lactate and creatinine.*
