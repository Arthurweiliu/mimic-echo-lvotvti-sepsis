source("config.R")
suppressMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(x))
d <- fread(file.path(DATA, "h1_data.csv"), colClasses=list(character=c('stay_id','subject_id')))
d[, icu_intime := as.POSIXct(icu_intime, tz='UTC')]
for (v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti','vti_diam','lvef'))
  d[, (v) := num(get(v))]
lk <- fread(file.path(DATA, "linked_cohort.csv"), colClasses=list(character='stay_id'))
lk <- lk[!duplicated(stay_id), .(stay_id, echo_dt_minus_intime_h)]
d <- merge(d, lk, by='stay_id', all.x=TRUE, sort=FALSE)
d[, echo_h := num(echo_dt_minus_intime_h)]
ic <- fread(file.path(DATA, "sepsis3_icustays.csv"), colClasses=list(character=c('subject_id','stay_id')))
ic[, intime := as.POSIXct(intime, tz='UTC')]; setorder(ic, subject_id, intime)
first_ids <- ic[!duplicated(subject_id), stay_id]
d[, start := echo_h/24]

n0 <- nrow(d); cat('first-ICU patients with Sepsis-3 + TTE linkage:', n0, '\n')
v <- !is.na(d$vti) & d$vti > 0
cat('with structured LVOT-VTI:', sum(v), '\n')
w <- v & !is.na(d$echo_h) & d$echo_h >= 0 & d$echo_h <= 72
cat('  index TTE 0-72 h:', sum(w), '\n')
a <- w & !is.na(d$fu_days) & d$fu_days > 0 & d$fu_days > d$start
cat('  alive at index TTE:', sum(a), '| died before TTE:', sum(w) - sum(a), '\n')
f <- a & d$stay_id %in% first_ids
cat('  true first ICU stay:', sum(f), '| deaths 90d:', sum(d$event90[f]), '\n')
cc <- f & !is.na(d$age) & !is.na(d$sex) & !is.na(d$sofa) & !is.na(d$charlson) &
  !is.na(d$lactate) & !is.na(d$creat) & d$sex %in% c('M','F')
cat('  complete covariates (blood-gas lactate):', sum(cc), '| deaths', sum(d$event90[cc]), '\n')
h1 <- f & !is.na(d$lvef) & d$lvef > 0
cat('  H1 (VTI+LVEF, no covariate requirement):', sum(h1), '| deaths', sum(d$event90[h1]), '\n')
h1cc <- h1 & !is.na(d$age) & !is.na(d$sofa) & !is.na(d$charlson) & !is.na(d$lactate) & !is.na(d$creat)
cat('  H1 complete covariates:', sum(h1cc), '| deaths', sum(d$event90[h1cc]), '\n')
h2 <- h1 & !is.na(d$vti_diam)
cat('  VTI+LVEF+LVOT diam (eCI computable), no covariate requirement:', sum(h2), '| deaths', sum(d$event90[h2]), '\n')
