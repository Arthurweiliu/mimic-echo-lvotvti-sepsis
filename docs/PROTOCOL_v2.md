# Protocol v2 — LVOT-VTI in sepsis: prognostic marker or incremental predictor?

**Version:** 2.0
**Date:** 2026-11-05
**Supersedes:** Protocol v1 and SAP amendment v1
**Status:** methods revision after independent audit; analysis to be re-run on the corrected cohort

---

## 1. Rationale

The first version linked every eligible ICU stay to any transthoracic echocardiogram (TTE) between 72 h before ICU admission and 72 h after ICU discharge. That window was wider than the stated “within 72 h of ICU admission”: 15.3% of the primary analysis set had the index TTE more than 72 h after admission, 7.7% after ICU discharge, and 60 patients died before the index TTE. Because follow-up began at ICU admission, the pre-TTE period was counted as if the exposure were known at baseline — a look-ahead (immortal-time) design.

This protocol fixes the design before re-running any analysis. The scientific question is unchanged: in adults with sepsis who survive to an early TTE, does LVOT-VTI carry prognostic information beyond conventional risk factors, and does it add incremental predictive value? The estimand and the time zero are now explicit.

---

## 2. Objectives and estimands

### 2.1 Primary objective
Estimate the association between LVOT-VTI measured at the index TTE and all-cause mortality within 90 days of ICU admission, among patients who survive to the TTE.

**Primary estimand.** In adults with Sepsis-3 in their first ICU stay, alive at an index TTE performed within 72 h after ICU admission and with LVOT-VTI measured, the hazard ratio of all-cause death by day 90 comparing LVOT-VTI < 17 cm with ≥ 17 cm, with follow-up beginning at the index TTE (left truncation) and adjustment for age, sex, SOFA, Charlson index, lactate and creatinine.

### 2.2 Co-primary objective
The same contrast in a **72-h landmark** analysis: patients alive at 72 h after ICU admission with an index TTE in [0, 72 h], followed from 72 h to death or day 90.

### 2.3 Secondary objectives
1. Continuous dose–response (per +5 cm; restricted cubic spline).
2. Comparison of per-beat (VTI, stroke volume) and per-minute (cardiac output, cardiac index) measures.
3. Incremental predictive value over the clinical model (ΔAUC, NRI, IDI, decision-curve net benefit).
4. Effect modification (atrial fibrillation, tachycardia, aortic valve disease, image quality), reported as exploratory with FDR control.
5. Calibration and internal validation of the clinical model with and without LVOT-VTI.

---

## 3. Design and setting

Retrospective cohort study using MIMIC-IV v3.1 and MIMIC-IV-ECHO v1.0.1, linked by `subject_id`/`hadm_id`/`stay_id`. The study is reported according to STROBE, with the RECORD extension for routinely collected data and TRIPOD for the prognostic-model components. The protocol is written as a target-trial emulation: eligibility at ICU admission, a fixed exposure window, an explicit time zero, and a defined contrast.

---

## 4. Eligibility

### 4.1 Inclusion
1. Age ≥ 18 years.
2. Sepsis-3 (suspected infection plus SOFA ≥ 2) at the index ICU admission.
3. **First ICU stay per patient**, ordered by `intime` across *all* of that patient's ICU stays in the source census, whether or not the stay contains an echocardiogram. A patient whose first ICU stay has no eligible TTE is excluded; a TTE in a later ICU stay is **not** admitted. "First" means the first ICU stay with Sepsis-3; the alternative definition (first ICU stay of any kind) is reported as a sensitivity analysis.
   *Implementation note:* the predicate is `stay_id == stay_id[which.min(intime)]` over all census rows for the subject. The v1 link retained the earliest stay **that has a TTE**, which admitted 346 of 3,919 patients (8.8%) in a later stay; that rule does not satisfy this criterion.
4. Index TTE of type “TTE” with a structured LVOT-VTI measurement, performed **between ICU admission and ICU admission + 72 h**.
5. Alive at the time of the index TTE.
6. Non-zero follow-up (`fu_days > start`, i.e. alive at the index TTE).
   The primary analysis does **not** require complete covariates. Eligibility is clinical; the primary population is "all eligible stays". Covariates are handled by multiple imputation (§9.8), and the complete-case analysis is reported as a sensitivity analysis with its own n, event count and flow line. The proportion excluded by the complete-case restriction and their mortality are reported, because that restriction is not ignorable: in the corrected cohort 26.1% of eligible patients lack a lactate and have lower mortality (25.7% vs 35.8%).

### 4.2 Index TTE rule
- The index TTE is the **earliest** TTE whose `measurement_datetime` lies in `[icu_intime, icu_intime + 72 h]` **and which carries a structured LVOT-VTI measurement**. A study without a VTI measurement is not an index study; a patient whose only in-window study lacks a VTI is excluded as "no LVOT-VTI", and the number of patients who had a later in-window study with a VTI is reported separately.
- If two studies share the same `measurement_datetime` to the second, the one with a VTI is retained; if both have one, the lower `measurement_id` is retained. The tie-break is deterministic and reproducible from the query.
- All timing is computed from `measurement_datetime` (the structured measurement timestamp); report dates are never used.
- TTEs before ICU admission and after 72 h are **not** eligible in the primary analysis.
- A TTE performed after death is a data error and is excluded; the count is reported.

### 4.3 Exclusion
- Sepsis not confirmed; age < 18; subsequent ICU stays; no TTE in the window; no LVOT-VTI; death before the index TTE; missing primary covariates; zero follow-up.
- Invalidating measurement conditions (cardiac tamponade, moderate/large pericardial effusion, LVOT obstruction, valve vegetation) are excluded from the primary analysis and examined in a sensitivity analysis, as in v1.

### 4.4 Cohort flow
A predicate-to-count table will be produced for every step:
`Sepsis-3 ICU stays → first ICU stay → TTE in [0,72 h] → LVOT-VTI available → complete covariates → Set A → Set C`.
Every count must be reproducible by a standalone query; the final set must close exactly.

---

## 5. Time zero and follow-up

### 5.1 Primary analysis
- **Time scale:** days from ICU admission.
- **Entry (left truncation):** the date/time of the index TTE (`start = (echo_time − ICU_intime)/24`).
- **Exit:** death, or day 90 after ICU admission, whichever first.
- **Event:** all-cause death by day 90.
- **Model:** Cox proportional-hazards model with `Surv(start, stop, event)`; patients contribute risk only from the index TTE onward. Pre-TTE person-time is not attributed to either exposure group.
- Patients with `start ≥ stop` (death before or at the TTE) are excluded and counted.

### 5.2 Co-primary landmark
- Restrict to patients alive at 72 h with `start ≤ 3` days.
- Time scale from ICU admission; entry at day 3 (`Surv(3, stop, event)`); exit at death or day 90.
- This is a **different estimand** (survivors to day 3; 239 deaths in the first 72 h are excluded by design), not a robustness check on the same population. Report n, deaths and mortality alongside the primary analysis.
- Both §5.1 and §5.2 remove immortal time, by different mechanisms (left truncation vs landmark). Agreement between them is evidence about the entry-time assumption, not about immortal time. If the two disagree materially, report both and do not select the more favourable one.

### 5.3 Sensitivity time-zero analyses
- Time zero at the index TTE; outcome = death within 90 days **of the TTE**.
- Time zero at ICU admission with a 7-day landmark.
- Exclusion of patients whose TTE occurred after 48 h.

---

## 6. Exposure and derived variables

- **Primary exposure:** LVOT-VTI (cm), measured at the index TTE.
- **Contrasts:** < 17 vs ≥ 17 cm (primary); < 17 vs > 20 cm; per +5 cm; quintiles; restricted cubic spline.
- **Derived:** SV = VTI × π × (LVOT diameter/2)²; CO = SV × echo-time heart rate / 1000; CI = CO / BSA.
- **Physiological ranges:** SV 5–250 mL; CO 1–15 L/min; CI 0.5–8 L/min/m².
- **Set C** requires a recorded LVOT diameter, an echo-time heart rate, and a BSA. The exclusion breakdown is reported (not only “missing diameter”).

### 6.1 Same-TTE co-measurement and the derived-CI caveat
- VTI, LVEF, SV, CO and CI are taken from the **same index TTE**; no patient contributes measures from different studies.
- CI is **echocardiography-derived** (`SV × HR / BSA`), not invasively measured. It is labelled **eCI** throughout and is never compared with invasively measured CI.
- Two analysis sets are defined:
  - **H1 (primary head-to-head):** same-TTE **VTI + LVEF** available.
  - **H2 (secondary):** same-TTE **VTI + LVEF + eCI** available.
  The availability matrix (VTI, LVEF, eCI) and the mortality of included versus excluded patients are reported.
- VTI and eCI are algebraically linked (`eCI = VTI × LVOT area × HR / BSA`). They are **not independent**; the primary comparison is VTI vs LVEF, and the joint VTI + eCI model is a sensitivity analysis with collinearity diagnostics (correlation matrix, VIF).

---

## 7. Outcome

- **Primary:** all-cause mortality within 90 days of ICU admission.
- **Ascertainment:** MIMIC-IV death registry; administrative censoring at day 90; post-discharge deaths retained.
- **Secondary:** 28-day mortality; in-hospital mortality; AKI (KDIGO ≥ 1), AKI stage ≥ 3, RRT (reported, not used for the primary claim).
- No cause-of-death analysis (not recorded).

---

## 8. Covariates

Pre-specified adjustment set: age, sex, first-day SOFA, Charlson comorbidity index, first-day lactate, first-day creatinine.

Not adjusted in the primary model (to avoid over-adjustment / mediation): vasopressor dose, ventilation, LVEF, fluid balance. These are reported descriptively and examined in sensitivity analyses.

---

## 9. Statistical analysis

### 9.1 Primary model
```
log h_i(t) = log h_0(t) + β·I(VTI_i < 17) + γ′Z_i,   t ≥ start_i
```
Cox partial likelihood with left truncation (`Surv(start, stop, event)`), Efron ties. Report HR, 95% CI, p; Kaplan–Meier from the TTE; proportional-hazards diagnostics for every covariate (Schoenfeld residuals, global test) and a time-split model at 7 days. **Pre-specified action when the exposure violates proportional hazards (Schoenfeld p < 0.05):** the single HR is reported as a *time-averaged* estimate over days 0–90 and is accompanied by episode-specific HRs for `(start, 7]` and `(7, 90]` with their confidence intervals; the abstract may quote a single HR only together with the statement that the association was concentrated in the first week. A single HR is never reported as if the hazard ratio were constant. (In the corrected cohort the Schoenfeld test for VTI < 17 gives p = 0.010; episode-specific HRs are 1.37 (1.15–1.63) and 1.14 (0.99–1.32).)

### 9.2 Dose–response
- Continuous per +5 cm.
- Restricted cubic spline of **log VTI** with 3 df (two internal knots at the 33rd and 67th percentiles), reference = cohort median.
- Likelihood-ratio test: spline vs linear (non-linearity) and spline vs covariate-only (overall).
- Plot HR (95% CI) vs VTI; mark the 5th/95th percentiles; do not describe the curve as monotonic unless it is monotonic over the full observed range; report any upturn in the top 5%.

### 9.3 Per-beat vs per-minute
- Same adjustment set; VTI, SV, CO, CI and HR as separate linear exposures.
- Log-transformed models for the algebraic identity: `log CO = log SV + log HR − log 1000`; show identical −2 log-likelihood and coefficients on the log scale.
- Do not claim exact equivalence for untransformed models; report Δ−2LL.
- Spline models for CO and CI with the same knots and reference; test non-linearity; sensitivity restricted to CO 3–10 L/min and CI 1.5–5 L/min/m².

### 9.4 Calibration
- **Predicted 90-day risk, conditional on entry.** Because follow-up is left-truncated at the index TTE, the model estimates the hazard only from `start_i` onward. The predicted risk must be conditional on survival to entry:
  `p_i = 1 − [S_0(90)/S_0(start_i)]^{exp(lp_i)} = 1 − exp(−[H_0(90) − H_0(start_i)]·exp(lp_i))`,
  with `S_0` (or `H_0`) evaluated at `t = 90` and `t = start_i` from the fitted delayed-entry model. The unconditional form `1 − S_0(90)^{exp(lp_i)}` is **not** valid here: it charges every patient the hazard accrued before their own entry and produces a spurious negative calibration intercept (in the corrected cohort: −0.089 instead of −0.002).
- **Required check:** before reporting slope, intercept or Brier score, verify that the mean predicted risk equals the observed 90-day mortality to within 0.002; report both numbers.
- Because a logistic recalibration model ignores censoring, report the number of patients censored alive before day 90. If that number is non-zero, logistic recalibration is replaced by an IPCW or Kaplan–Meier-based calibration estimate.
- Calibration slope and intercept from logistic recalibration (`event ~ qlogis(p)`); ideal 1 and 0.
- Deciles of predicted risk: observed vs predicted with 95% binomial CIs; state the number of deciles outside their CI.
- Brier score.
- Bootstrap optimism-corrected calibration slope (300 resamples), with the spline knots fixed at the full-cohort 33rd/67th percentiles in every resample.

### 9.5 Discrimination and internal validation
- Harrell C-index; bootstrap optimism-corrected C-index (400 resamples).
- ΔC-index between the clinical model and the clinical model + LVOT-VTI; bootstrap percentile CI of the apparent ΔC and the optimism-corrected ΔC, clearly labelled.
- No external validation is possible; state this as a limitation.

### 9.6 Incremental value
- Logistic models for 28-day and 90-day mortality.
- ΔAUC; category-free NRI; risk-category NRI (<10%, 10–25%, ≥25%); IDI; decision-curve net benefit.
- Pre-defined meaningful increment: ΔAUC ≥ 0.01, or risk-category NRI ≥ 0.10 with a CI excluding 0, or a net-benefit advantage with a CI excluding 0 at a clinically relevant threshold. Otherwise “no meaningful incremental value”.
- Every metric is computed in the **primary cohort** (n and event count as in §5.1), not in a wider complete-case set; the analysis set, event count and number of observations are stated with each metric.
- **Head-to-head models** are fitted in the same n: clinical; +VTI; +LVEF; +eCI; +VTI+LVEF; +VTI+eCI; +LVEF+eCI; all. Compare AIC, −2 log-likelihood, C-index, calibration, NRI/IDI and decision-curve net benefit; use paired bootstrap for ΔC. VTI, LVEF and eCI are **not imputed**; covariates are handled by MI (§9.8).

### 9.7 Effect modification
- Pre-specified: atrial fibrillation/flutter, tachycardia (echo HR > 100), aortic stenosis ≥ moderate, image quality.
- Additional: aortic/mitral/tricuspid regurgitation, ICU mean HR.
- One family, Benjamini–Hochberg FDR; report q values to 3 decimals; all exploratory.

### 9.8 Missing data
- Multiple imputation by chained equations (m = 20) for the six covariates, using all analysis variables plus auxiliary variables (SOFA components, vasopressor use, ventilation, ICU length of stay) in the imputation model; Rubin’s rules.
- **Multiple imputation is the primary analysis; complete-case analysis is the sensitivity analysis.** The two answer different questions and both are reported with their own cohort-flow line.
- The exposure (LVOT-VTI) is **never** imputed. The estimand is defined “with LVOT-VTI measured”; imputing the exposure would change the population rather than recover it (40.2% of the linked cohort has no VTI).
- The outcome is **never** imputed: it is registry-ascertained and complete to day 90, and imputing it would add assumptions without adding information.
- Report the missingness mechanism and the proportion missing for each variable.

### 9.9 Sensitivity analyses
1. Window 0–24 h; 0–48 h; [−24 h, +72 h]. Each restricts the population and is reported with its own n and event count; it is not a sensitivity analysis of the same estimand.
2. True first ICU stay (all census rows, not “first stay with a TTE”); report the 346 excluded in the current pipeline.
3. Exclude death before TTE (already in primary; report the count, 53).
4. 72-h landmark (co-primary).
5. Time zero at TTE; outcome 90 days from TTE.
6. Exclude invalidating measurement conditions; exclude moderate-or-worse aortic valve disease; adjust for atrial arrhythmia, valve disease and image quality.
7. In-hospital death only (competing-risk estimator specified; discharge alive is a competing event).
8. Inverse-probability-of-selection weighting for TTE vs no TTE, with a stated target population and positivity check.
9. All eligible stays with patient-clustered standard errors — presented as a different population, not an SE variant.
10. Measurement-error simulation (5% error in VTI and LVOT diameter).
11. E-value for unmeasured confounding, explicitly distinguished from selection on echo availability.
12. Exclude the top/bottom 1% of VTI.
13. Add LVEF to the clinical model, so the incremental value of LVOT-VTI is not measured against an artificially weak benchmark.

### 9.10 Multiplicity
- Primary and co-primary: two tests; report both; no alpha adjustment (pre-specified dual primary).
- Secondary and exploratory: no confirmatory claims; FDR for the effect-modification family only.
- All secondary results labelled exploratory.

### 9.11 Software and reproducibility
R 4.5.2 (`survival`, `splines`); code in a public repository; a frozen analysis dataset and a hash of the manuscript; all scripts run from raw CSVs; no manual edits to numbers.

---

## 10. Bias mitigation

| Bias | Design/analysis safeguard |
|---|---|
| Immortal time / look-ahead | Fixed [0, 72 h] window; left truncation at TTE; 72-h landmark; exclude death before TTE; report the count |
| Confounding by indication (TTE) | Adjustment; IPW for selection into the TTE cohort; report SMDs; no causal claim |
| Selection (first ICU stay) | One record per patient; report subsequent-stay exclusions; cluster-robust sensitivity |
| Measurement error | Single TTE acknowledged; 5% error simulation; RCS to avoid cut-point dependence |
| Missing data | MI (m = 20) primary; complete-case sensitivity; exposure and outcome never imputed; report |
| Outcome ascertainment | Registry + in-hospital death; restrict to in-hospital death as sensitivity |
| Multiplicity | Pre-specified primary/co-primary; FDR for modifiers; exploratory labels |
| Functional-form artefact | RCS on log exposure; report linear and spline; no “monotonic” claim without full-range check |

---

## 11. Reporting and audit

- STROBE checklist; RECORD items for database studies; TRIPOD for the prognostic-model section.
- A **predicate-to-count** table for every cohort step, produced before submission.
- An **independent audit** of the frozen manuscript and code: re-derive the cohort, the primary model, the spline p-values and the calibration statistics from the raw CSVs; no shared code with the primary analyst.
- Freeze the manuscript (SHA-256) before audit; any change after the freeze creates a new version and a new audit.
- The audit report is an internal document, not a journal supplement.

---

## 12. Changes from v1 (amendment log)

| # | v1 | v2 | Reason |
|---|---|---|---|
| 1 | TTE window `[intime−72 h, outtime+72 h]` | `[intime, intime+72 h]` | The v1 window was wider than the stated criterion; it admitted post-discharge echoes and 60 deaths before the echo |
| 2 | Time zero at ICU admission, exposure at any time in the window | Time zero at ICU admission; **left truncation at the index TTE** | Removes immortal time without discarding the ICU-admission time scale |
| 3 | No co-primary landmark | 72-h landmark co-primary | Pre-specified design that removes immortal time by construction |
| 4 | Linear models only for CO/CI | RCS on log exposure for all flow measures | The linear CO null was a functional-form artefact |
| 5 | Calibration not reported | Slope, intercept, deciles with CIs, Brier | Editor and TRIPOD requirements |
| 6 | Internal validation not reported | Bootstrap optimism-corrected C-index and ΔC | Editor requirement |
| 7 | “Monotonic per-beat” | “Inverse across the clinical range; check the top 5%” | The SV curve rises at the top 5% |
| 8 | “Equivalent” CO/SV models | Exact only on the log scale; untransformed Δ−2LL reported | Arithmetic and scale correctness |
| 9 | 600 excluded from Set C attributed to missing diameter | Full breakdown (diameter, BSA, HR, range) | Only 298/600 lacked a diameter |
| 10 | No cohort audit | Predicate-to-count table + independent audit | Prevent recurrence |

---

## Appendix A. Current corrected estimates (planning only, not for publication)

| Analysis | n | Deaths | Estimate |
|---|---:|---:|---|
| Primary, left-truncated Cox, VTI < 17 vs ≥ 17 (all in-window, complete covariates) | 3,919 | 1,403 | HR 1.23 (1.10–1.37) |
| **Protocol-compliant true first ICU stay** | **3,573** | **1,262** | **HR 1.26 (1.12–1.42)** |
| < 17 vs > 20 cm | 3,919 | 1,403 | 1.25 (1.10–1.41) |
| per +5 cm | 3,919 | 1,403 | 0.957 (0.919–0.996) |
| 72-h landmark | 3,680 | 1,164 | 1.20 (1.06–1.36) |
| Time-split, days (start, 7] | 3,919 | — | 1.37 (1.15–1.63) |
| Time-split, days (7, 90] | 3,919 | — | 1.14 (0.99–1.32) |
| Set C, VTI +5 cm | 3,453 | 1,212 | 0.944 (0.903–0.987) |
| Set C, SV +10 mL | 3,453 | 1,212 | 0.972 (0.947–0.997) |
| Set C, CO +1 L/min | 3,453 | 1,212 | 1.003 (0.975–1.031) |
| Set C, CI +1 | 3,453 | 1,212 | 1.046 (0.993–1.102) |
| Set C, HR +10 bpm | 3,453 | 1,212 | 1.062 (1.031–1.093) |
| Spline non-linearity, VTI | — | — | p = 0.0012 |
| Spline non-linearity, SV | — | — | p = 0.00012 |
| Spline non-linearity, CO | — | — | p = 0.0044 |
| Spline non-linearity, CI | — | — | p = 0.0026 |
| Calibration, entry-conditional (Set A) | 3,919 | 1,403 | slope 0.996; intercept −0.002; mean predicted 0.3579 vs observed 0.3580 |
| Incremental value (Set A) | 3,919 | 1,403 | ΔAUC +0.0012; ΔC +0.001 |

These numbers are provisional and must be reproduced by the audit before they enter a manuscript. The 3,919 row includes 346 patients (8.8%) whose index stay is not their first ICU stay; the protocol-compliant row restricts to true first stays and is the estimate that should be quoted.
