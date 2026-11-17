### Table S3. Linear versus non-linear functional form for every same-study measure (multiple imputation)

| Measure | Contrast | n | Deaths | Linear HR (95% CI) | p linear | Non-linearity p | Overall p | delta AIC (spline - linear) |
|---|---|---:|---:|---|---:|---:|---:|---:|
| LVOT-VTI | per 5 | 4,851 | 1,608 | 0.93 (0.90-0.97) | <0.001 | 0.00093 | 3.1e-08 | +10.0 |
| Stroke volume | per 10 | 4,579 | 1,508 | 0.97 (0.95-0.99) | 0.003 | 0.00019 | 2.3e-07 | +13.1 |
| Cardiac output | per 1 | 4,534 | 1,500 | 1.01 (0.99-1.03) | 0.440 | 0.0009 | 0.0024 | +10.0 |
| Cardiac index | per 1 | 4,307 | 1,413 | 1.06 (1.01-1.11) | 0.018 | 0.00039 | 0.00069 | +11.7 |
| Heart rate | per 10 | 4,803 | 1,598 | 1.08 (1.05-1.11) | <0.001 | 0.99 | 1.9e-07 | -4.0 |
| LVEF | per 10 | 4,369 | 1,436 | 0.97 (0.94-1.00) | 0.049 | 6.1e-12 | 5.4e-12 | +47.6 |

*Covariates imputed (m = 20) within each measure-specific cohort; the markers themselves are not imputed. The linear coefficient is fitted on the untransformed exposure (the quantity under test) and pooled by Rubin's rules; the non-linearity test compares a 3-df natural spline with a linear term on the same scale (log for the flow measures, untransformed for LVEF and heart rate) and is reported as the median p across the 20 imputations. A positive delta AIC means the spline fits better.*
