## 106_strict_first_icu.R — (1) how many cohort patients had an earlier ICU stay of any kind,
## (2) sensitivity analysis restricted to the first ICU stay of any kind,
## (3) verification of the "death before the index study" records against admissions.deathtime
source("config.R")
suppressMessages(library(data.table)); suppressMessages(library(survival)); suppressMessages(library(mice))
num <- function(x) suppressWarnings(as.numeric(x))
set.seed(20261105)

## ---------- icustays: rank of every ICU stay within its patient ----------
icu <- fread(file.path(DATA, "icustays.csv"),
             colClasses=list(character=c('subject_id','hadm_id','stay_id')))
icu[, intime := as.POSIXct(intime, format='%Y/%m/%d %H:%M:%S', tz='UTC')]
setorder(icu, subject_id, intime)
icu[, rank_any := seq_len(.N), by = subject_id]
cat('ICU stays:', nrow(icu), '| patients:', uniqueN(icu$subject_id), '\n')

## ---------- rebuild the 4,816 cohort exactly as in 66_counts.R ----------
d <- fread(file.path(DATA, "h1_data.csv"), colClasses=list(character=c('stay_id','subject_id')))
for (v in c('fu_days','event90','vti','lvef','vti_diam')) d[, (v) := num(get(v))]
lk <- fread(file.path(DATA, "linked_cohort.csv"), colClasses=list(character='stay_id'))
lk <- lk[!duplicated(stay_id), .(stay_id, echo_dt_minus_intime_h)]
d <- merge(d, lk, by='stay_id', all.x=TRUE, sort=FALSE)
d[, echo_h := num(echo_dt_minus_intime_h)]
s3 <- fread(file.path(DATA, "sepsis3_icustays.csv"),
            colClasses=list(character=c('subject_id','stay_id')))
s3[, intime := as.POSIXct(intime, tz='UTC')]; setorder(s3, subject_id, intime)
first_s3 <- s3[!duplicated(subject_id), stay_id]
d[, start := echo_h / 24]
A <- d[!is.na(vti) & vti > 0 & !is.na(fu_days) & fu_days > 0 & !is.na(echo_h) &
       echo_h >= 0 & echo_h <= 72 & fu_days > start & stay_id %in% first_s3]
cat('analytic cohort:', nrow(A), 'deaths', sum(A$event90), '\n')

A <- merge(A, icu[, .(stay_id, rank_any, intime_icu = intime)], by='stay_id', all.x=TRUE)
cat('cohort rows with an icustays match:', sum(!is.na(A$rank_any)), '\n')

## ---------- (1) how many entered through a later ICU stay ----------
tab <- A[, .N, by = .(later = rank_any > 1)]
print(tab)
cat(sprintf('\nentered through a later ICU stay (any cause): %d of %d (%.1f%%)\n',
            sum(A$rank_any > 1), nrow(A), 100 * mean(A$rank_any > 1)))
cat('  their 90-day mortality: ', sprintf('%.1f%%', 100 * mean(A$event90[A$rank_any > 1])),
    '| first-stay patients: ', sprintf('%.1f%%', 100 * mean(A$event90[A$rank_any == 1])), '\n', sep='')

## ---------- (2) sensitivity: first ICU stay of any kind ----------
A1 <- A[rank_any == 1]
cat(sprintf('\nstrict first-ICU-stay cohort: n = %d, deaths = %d\n', nrow(A1), sum(A1$event90)))

A1 <- merge(A1, fread('lactate_v2.csv', colClasses=list(character=c('stay_id','subject_id')))[
  , .(stay_id, lactate_filled2)], by='stay_id', all.x=TRUE, sort=FALSE)
A1[, vlow := as.numeric(vti < 17)]
AUX <- c('sapsii','mbp_min','hr_mean','temp_mean','spo2_min','bun','hgb','plt','wbc','bili',
         'inr','norepi','ventN','rrtN','vti_diam')
tf <- function(x) { s <- trimws(as.character(x)); s[is.na(s)] <- ''; ifelse(grepl('^t', tolower(s)), 1, 0) }
A1[, `:=`(ventN = tf(vent), rrtN = tf(rrt))]
COV <- c('age','sex','sofa','charlson','lactate_filled2','creat')
X <- as.data.frame(A1[, c('start','fu_days','event90','vlow', COV, AUX), with=FALSE])
meth <- rep('', ncol(X)); names(meth) <- names(X)
for (v in c(COV, AUX)) if (sum(is.na(X[[v]])) > 0) meth[v] <- 'pmm'
sf <- survfit(Surv(start, fu_days, event90) ~ 1, data=X)
X$nag <- approx(sf$time, sf$cumhaz, X$fu_days, rule=2, ties='ordered')$y
X$nag[is.na(X$nag)] <- 0
meth <- c(meth, nag='')
imp <- mice(X, m=20, maxit=10, method=meth, printFlag=FALSE, seed=20261105)
fits <- lapply(seq_len(imp$m), function(k)
  coxph(Surv(start,fu_days,event90) ~ vlow + age + sex + sofa + charlson + lactate_filled2 + creat,
        data=complete(imp,k)))
pl <- summary(pool(as.mira(fits)), conf.int=TRUE, exponentiate=TRUE)
b <- pl[grep('^vlow', pl$term), ]
cat(sprintf('SENSITIVITY (first ICU stay of any kind): HR %.3f (%.3f-%.3f)\n',
            b$estimate, b[["2.5 %"]], b[["97.5 %"]]))
write.csv(data.frame(analysis='First ICU stay of any kind (sensitivity)', n=nrow(A1),
                     deaths=sum(A1$event90), hr=b$estimate, lo=b[["2.5 %"]], hi=b[["97.5 %"]]),
          'TableS8_strict_first_icu.csv', row.names=FALSE)

## (3) death timing is checked in 106b_death_timing.R
