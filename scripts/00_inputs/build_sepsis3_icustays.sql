-- build_sepsis3_icustays.sql — export Sepsis-3 ICU stays (the census that seeds the H1 cohort).
-- Same join as inside build_cohort_covariates.sql; produced 32,971 rows in the published analysis.
-- Run psql FROM YOUR DATA DIRECTORY: the \copy output is relative to the working directory.
--   cd /path/to/data && psql -d mimic -f build_sepsis3_icustays.sql
\copy (select d.subject_id, d.stay_id, d.hadm_id, d.icu_intime as intime, d.icu_outtime as outtime, d.los_icu as los from mimiciv_derived.icustay_detail d join mimiciv_derived.sepsis3 s on s.stay_id = d.stay_id and s.sepsis3) to 'sepsis3_icustays.csv' csv header
