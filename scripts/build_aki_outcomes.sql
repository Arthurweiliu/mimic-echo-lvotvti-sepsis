-- build_aki_outcomes.sql — export AKI / RRT outcomes from PostgreSQL.
-- Run psql FROM YOUR DATA DIRECTORY: the \copy output is relative to the working directory.
--   cd /path/to/data && psql -d mimic -f build_aki_outcomes.sql
\copy (select d.stay_id, max(ks.aki_stage) as max_aki_stage_7d, max(ks.aki_stage_crrt) as max_aki_crrt_7d, bool_or(ks.aki_stage >= 1) as aki_any_7d, bool_or(ks.aki_stage >= 3) as aki3_7d, bool_or(ks.aki_stage_crrt >= 1) as crrt_7d from mimiciv_derived.icustay_detail d join mimiciv_derived.sepsis3 s on s.stay_id = d.stay_id and s.sepsis3 join mimiciv_derived.kdigo_stages ks on ks.stay_id = d.stay_id and ks.charttime between d.icu_intime and d.icu_intime + interval '7 day' group by d.stay_id) to 'aki_outcomes.csv' csv header
