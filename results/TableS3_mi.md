### Table S3. Linear versus non-linear functional form for every same-study measure (multiple imputation)

| Measure | Contrast | n | Deaths | Linear HR (95% CI) | p linear | Non-linearity p | Overall p | delta AIC (spline - linear) |
|---|---|---:|---:|---|---:|---:|---:|---:|
| LVOT-VTI | per 5 | 4,816 | 1,573 | 0.93 (0.90-0.97) | <0.001 | 0.00085 | 3.8e-08 | +10.1 |
| Stroke volume | per 10 | 4,547 | 1,476 | 0.97 (0.95-0.99) | 0.006 | 0.00019 | 5.8e-07 | +13.2 |
| Cardiac output | per 1 | 4,502 | 1,468 | 1.01 (0.99-1.04) | 0.302 | 0.00057 | 0.0018 | +10.9 |
| Cardiac index | per 1 | 4,276 | 1,382 | 1.06 (1.02-1.11) | 0.008 | 0.00027 | 0.00037 | +12.5 |
| Heart rate | per 10 | 4,768 | 1,563 | 1.08 (1.05-1.11) | <0.001 | 0.97 | 2.1e-07 | -3.9 |
| LVEF | per 10 | 4,336 | 1,403 | 0.97 (0.94-1.00) | 0.065 | 1.9e-12 | 2.1e-12 | +50.0 |

*Covariates imputed (m = 20) within each measure-specific cohort; the markers themselves are not imputed. The linear coefficient is fitted on the untransformed exposure (the quantity under test) and pooled by Rubin's rules; the non-linearity test compares a 3-df natural spline with a linear term on the same scale (log for the flow measures, untransformed for LVEF and heart rate) and is reported as the median p across the 20 imputations. A positive delta AIC means the spline fits better.*
