## 108_build_precise_outcomes.R — rebuild fu_days / event90 / death28 from precise timestamps
## deathtime (admissions) where the patient died in hospital; otherwise the recorded date of death.
setwd('F:/openclaw_workspace/H1_revision')
suppressMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(x))

d <- fread('F:/openclaw_workspace/sepsis_echo_lit/h1_data.csv', colClasses=list(character=c('stay_id','subject_id')))
s3 <- fread('F:/openclaw_workspace/sepsis_echo_lit/sepsis3_icustays.csv',
            colClasses=list(character=c('subject_id','stay_id','hadm_id')))
s3[, intime := as.POSIXct(intime, tz='UTC')]
adm <- fread('F:/openclaw_workspace/sepsis_echo_lit/admissions.csv',
             colClasses=list(character=c('subject_id','hadm_id')),
             select=c('subject_id','hadm_id','admittime','dischtime','deathtime','hospital_expire_flag'))
adm[, deathtime := as.POSIXct(deathtime, format='%Y/%m/%d %H:%M:%S', tz='UTC')]
adm[, dischtime := as.POSIXct(dischtime, format='%Y/%m/%d %H:%M:%S', tz='UTC')]

d <- merge(d, s3[, .(stay_id, hadm_id, intime)], by='stay_id', all.x=TRUE)
d <- merge(d, adm[, .(subject_id, hadm_id, deathtime, dischtime, hospital_expire_flag)],
           by=c('subject_id','hadm_id'), all.x=TRUE)

cov <- fread('F:/openclaw_workspace/sepsis_echo_lit/cohort_covariates.csv',
             colClasses = list(character = c('subject_id','hadm_id','stay_id')),
             select = c('stay_id','dod'))
cov <- cov[!duplicated(stay_id)]
d <- merge(d, cov, by = 'stay_id', all.x = TRUE)
d[, dod_ts := as.POSIXct(as.Date(dod), tz='UTC')]
## death timestamp: precise deathtime first, then in-hospital dischtime, then the recorded date of death
d[, death_ts := fifelse(!is.na(deathtime), deathtime,
                 fifelse(hospital_expire_flag == 1 & !is.na(dischtime), dischtime, dod_ts))]
d[, death_src := fifelse(!is.na(deathtime), 'deathtime',
                  fifelse(hospital_expire_flag == 1 & !is.na(dischtime), 'dischtime', 'dod'))]

d[, fu_precise := as.numeric(difftime(death_ts, intime, units='days'))]
d[, event90_v2 := as.numeric(!is.na(fu_precise) & fu_precise <= 90)]
d[, fu_days_v2 := pmin(fu_precise, 90)]
d[is.na(fu_days_v2) | fu_days_v2 <= 0, `:=`(fu_days_v2 = 90, event90_v2 = 0)]
d[, death28_v2 := as.numeric(!is.na(fu_precise) & fu_precise <= 28)]
cat('deaths by source of timestamp:\n'); print(d[event90_v2 == 1, .N, by = death_src])
cat('\n90-day events: old', sum(d$event90, na.rm=TRUE), '| new', sum(d$event90_v2), '\n')
cat('rows where the event changed:', d[, sum(event90 != event90_v2, na.rm=TRUE)], '\n')
cat('max |fu_days new - old| among deaths:', round(max(abs(d$fu_days_v2 - d$fu_days)[d$event90_v2 == 1]), 2), 'days\n')

d[, `:=`(fu_days = fu_days_v2, event90 = event90_v2, death28 = death28_v2)]
out <- d[, c(names(fread('F:/openclaw_workspace/sepsis_echo_lit/h1_data.csv', nrows=0)), 'death_src'), with=FALSE]
keep <- names(out)[names(out) %in% names(d)]
fwrite(d[, ..keep], 'h1_data_v2.csv')
cat('\nwritten h1_data_v2.csv with', nrow(d), 'rows and', length(keep), 'columns\n')
