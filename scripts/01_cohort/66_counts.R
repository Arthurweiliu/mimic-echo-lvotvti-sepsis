## 66_counts.R — cohort and missingness counts quoted in the manuscript (no models)
source("config.R")
suppressMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(x))
d <- fread(file.path(DATA, 'h1_data.csv'), colClasses=list(character=c('stay_id','subject_id')))
for (v in c('fu_days','event90','vti','lvef','vti_diam','lactate')) d[, (v) := num(get(v))]
lk <- fread(file.path(DATA, "linked_cohort.csv"), colClasses=list(character='stay_id'))
lk <- lk[!duplicated(stay_id), .(stay_id, echo_dt_minus_intime_h)]
d <- merge(d, lk, by='stay_id', all.x=TRUE, sort=FALSE)
d[, echo_h := num(echo_dt_minus_intime_h)]
ic <- fread(file.path(DATA, "sepsis3_icustays.csv"), colClasses=list(character=c('subject_id','stay_id')))
ic[, intime := as.POSIXct(intime, tz='UTC')]; setorder(ic, subject_id, intime)
first_ids <- ic[!duplicated(subject_id), stay_id]
d[, start := echo_h/24]
inwin  <- !is.na(d$vti) & d$vti > 0 & !is.na(d$echo_h) & d$echo_h >= 0 & d$echo_h <= 72
alive  <- inwin & !is.na(d$fu_days) & d$fu_days > d$start
first  <- alive & d$stay_id %in% first_ids
A <- merge(d[first, ], fread('lactate_v2.csv', colClasses=list(character=c('stay_id','subject_id'))),
           by='stay_id', all.x=TRUE, suffixes=c('', '.v2'))

cat('first-ICU patients (Sepsis-3 + TTE):', nrow(d), '\n')
cat('with LVOT-VTI:', sum(!is.na(d$vti) & d$vti > 0), '\n')
cat('index TTE 0-72 h:', sum(inwin), '\n')
cat('alive at index TTE:', sum(alive), '| died before TTE:', sum(inwin) - sum(alive), '\n')
cat('true first ICU stay (PRIMARY COHORT):', nrow(A), '| deaths:', sum(A$event90),
    sprintf('| 90-day mortality %.1f%%', 100*mean(A$event90)), '\n')
cat('not first ICU stay (excluded from primary):', sum(alive) - nrow(A), '\n')
cat('blood-gas lactate present:', sum(!is.na(A$lactate)),
    '| missing:', sum(is.na(A$lactate)), sprintf('(%.1f%%)', 100*mean(is.na(A$lactate))), '\n')
cat('after chemistry supplement: present', sum(!is.na(A$lactate_filled2)),
    '| missing', sum(is.na(A$lactate_filled2)), sprintf('(%.1f%%)', 100*mean(is.na(A$lactate_filled2))), '\n')
cat('after nearest rule: missing', sum(is.na(A$lactate_nearest)), '\n')
cat('mortality by lactate availability (after supplement):\n')
print(A[, .(n = .N, deaths = sum(event90), mort = sprintf('%.1f%%', 100*mean(event90))),
        by = .(missing = is.na(lactate_filled2))])
cat('H1 (VTI + LVEF):', sum(!is.na(A$lvef) & A$lvef > 0), '| deaths',
    sum(A$event90[!is.na(A$lvef) & A$lvef > 0]), '\n')
cat('H2 (VTI + LVEF + LVOT diameter):', sum(!is.na(A$lvef) & A$lvef > 0 & !is.na(A$vti_diam)),
    '| deaths', sum(A$event90[!is.na(A$lvef) & A$lvef > 0 & !is.na(A$vti_diam)]), '\n')
h1 <- A[!is.na(lvef) & lvef > 0]; h2 <- A[!is.na(lvef) & lvef > 0 & !is.na(vti_diam)]
cat(sprintf('H1 mortality %.1f%%; without same-study LVEF n=%d (%.1f%%)\n',
            100*mean(h1$event90), nrow(A) - nrow(h1), 100*mean(A$event90[!(A$stay_id %in% h1$stay_id)])))
cat(sprintf('H2 mortality %.1f%%; not in H2 n=%d (%.1f%%)\n',
            100*mean(h2$event90), nrow(A) - nrow(h2), 100*mean(A$event90[!(A$stay_id %in% h2$stay_id)])))
cat('deaths by day 7:', sum(A$event90 == 1 & A$fu_days <= 7),
    '| deaths after day 7:', sum(A$event90 == 1 & A$fu_days > 7),
    '| alive at day 7:', sum(A$fu_days > 7), '\n')
cat('deaths by day 3 (excluded from the 72-h landmark):', sum(A$event90 == 1 & A$fu_days <= 3), '\n')
cat('VTI <17 or >20:', sum(A$vti < 17 | A$vti > 20), '| deaths', sum(A$event90[A$vti < 17 | A$vti > 20]), '\n')
