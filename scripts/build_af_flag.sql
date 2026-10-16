-- build_af_flag.sql — ICU-documented atrial fibrillation/ flutter from chartevents (Heart Rhythm)
-- plus ICD-coded AF/flutter per hospital admission.  Run from the directory holding _coh_stayids.csv.
--   psql -d mimiciv -f _build_af_flag.sql
\set ON_ERROR_STOP on
create temp table cs(stay_id bigint);
\copy cs from '_coh_stayids.csv' csv header

-- values of the Heart Rhythm item, for reference
\echo '=== Heart Rhythm value distribution (cohort) ==='
select lower(trim(value)) as rhythm, count(*) as n
from cs c join mimiciv_icu.chartevents ce on ce.stay_id = c.stay_id and ce.itemid = 220048
group by 1 order by 2 desc limit 25;

\echo '=== AF flag per stay (ICU chart) ==='
\copy (select c.stay_id, max(case when lower(ce.value) like '%fib%' or lower(ce.value) like '%afib%' or lower(trim(ce.value)) in ('af','a fib','afib') then 1 else 0 end) as af_icu from cs c join mimiciv_icu.chartevents ce on ce.stay_id = c.stay_id and ce.itemid = 220048 group by c.stay_id) to 'af_icu.csv' csv header

\echo '=== AF flag per hospital admission (ICD) ==='
\copy (select c.stay_id, max(case when d.icd_code like 'I48%' or d.icd_code like '427.3%' then 1 else 0 end) as af_icd from cs c join mimiciv_icu.icustays ie on ie.stay_id = c.stay_id join mimiciv_hosp.diagnoses_icd d on d.hadm_id = ie.hadm_id group by c.stay_id) to 'af_icd.csv' csv header
