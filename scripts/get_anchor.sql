-- get_anchor.sql — export anchor year groups from PostgreSQL.
-- Run psql FROM YOUR DATA DIRECTORY: the \copy output is relative to the working directory.
--   cd /path/to/data && psql -d mimic -f get_anchor.sql
\copy (select subject_id, anchor_year_group from mimiciv_hosp.patients) to 'anchor_year_group.csv' csv header
