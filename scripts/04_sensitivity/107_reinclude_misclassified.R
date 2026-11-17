## 107_reinclude_misclassified.R — patients excluded by the "alive at the index study" filter
## only because death is recorded to the day; re-included with precise deathtime.
source("config.R")
suppressMessages(library(data.table)); suppressMessages(library(survival)); suppressMessages(library(mice))
num <- function(x) suppressWarnings(as.numeric(x))
set.seed(20261105)

d <- fread(file.path(DATA, "h1_data.csv"), colClasses=list(character=c('stay_id','subject_id')))
for (v in c('fu_days','event90','vti','lvef','vti_diam','age','sofa','charlson','creat')) d[, (v) := num(get(v))]
lk <- fread(file.path(DATA, "linked_cohort.csv"), colClasses=list(character='stay_id'))
lk <- lk[!duplicated(stay_id), .(stay_id, echo_h = num(echo_dt_minus_intime_h))]
d <- merge(d, lk, by='stay_id', all.x=TRUE, sort=FALSE)
d[, start := echo_h / 24]

s3 <- fread(file.path(DATA, "sepsis3_icustays.csv"),
            colClasses=list(character=c('subject_id','stay_id','hadm_id')))
s3[, intime := as.POSIXct(intime, tz='UTC')]; setorder(s3, subject_id, intime)
first_s3 <- s3[!duplicated(subject_id), stay_id]
d <- merge(d, s3[, .(stay_id, hadm_id, intime)], by='stay_id', all.x=TRUE)

adm <- fread(file.path(DATA, "admissions.csv"),
             colClasses=list(character=c('subject_id','hadm_id')),
             select=c('subject_id','hadm_id','deathtime'))
adm <- adm[!is.na(deathtime) & deathtime != '']
adm[, deathtime := as.POSIXct(deathtime, format='%Y/%m/%d %H:%M:%S', tz='UTC')]
d <- merge(d, adm, by=c('subject_id','hadm_id'), all.x=TRUE)

## base eligibility: VTI, window, first Sepsis-3 stay
base <- d[!is.na(vti) & vti > 0 & !is.na(fu_days) & fu_days > 0 & !is.na(echo_h) &
          echo_h >= 0 & echo_h <= 72 & stay_id %in% first_s3]
cat('base (VTI, 0-72 h, first Sepsis-3 stay):', nrow(base), '\n')

## misclassified: fu_days <= start but the precise deathtime is after the study
mis <- base[fu_days <= start & !is.na(deathtime) & deathtime > intime + echo_h * 3600]
mis[, fu_precise := as.numeric(difftime(deathtime, intime, units='days'))]
cat('reclassified as alive at the study (precise deathtime):', nrow(mis), '\n')
cat(sprintf('  their corrected follow-up: median %.1f days (IQR %.1f to %.1f)\n',
            median(mis$fu_precise), quantile(mis$fu_precise, .25), quantile(mis$fu_precise, .75)))
cat('  deaths within 90 days:', sum(mis$fu_precise <= 90), '\n')

## corrected cohort = primary + reclassified
keep <- base[fu_days > start]
add <- copy(mis)
add[, `:=`(fu_days = pmin(fu_precise, 90), event90 = as.numeric(fu_precise <= 90))]
A2 <- rbind(keep[, .(stay_id, subject_id, fu_days, event90, vti, start, age, sex, sofa, charlson, creat,
                     sapsii, mbp_min, hr_mean, temp_mean, spo2_min, bun, hgb, plt, wbc, bili, inr, norepi,
                     vent, rrt, vti_diam)],
            add[, .(stay_id, subject_id, fu_days, event90, vti, start, age, sex, sofa, charlson, creat,
                    sapsii, mbp_min, hr_mean, temp_mean, spo2_min, bun, hgb, plt, wbc, bili, inr, norepi,
                    vent, rrt, vti_diam)])
cat(sprintf('\ncorrected cohort: n = %d, deaths = %d (primary 4,816 / 1,573)\n', nrow(A2), sum(A2$event90)))

A2 <- merge(A2, fread('lactate_v2.csv', colClasses=list(character=c('stay_id','subject_id')))[
  , .(stay_id, lactate_filled2)], by='stay_id', all.x=TRUE, sort=FALSE)
A2[, vlow := as.numeric(vti < 17)]
tf <- function(x) { s <- trimws(as.character(x)); s[is.na(s)] <- ''; ifelse(grepl('^t', tolower(s)), 1, 0) }
A2[, `:=`(ventN = tf(vent), rrtN = tf(rrt))]
AUX <- c('sapsii','mbp_min','hr_mean','temp_mean','spo2_min','bun','hgb','plt','wbc','bili',
         'inr','norepi','ventN','rrtN','vti_diam')
COV <- c('age','sex','sofa','charlson','lactate_filled2','creat')
X <- as.data.frame(A2[, c('start','fu_days','event90','vlow', COV, AUX), with=FALSE])
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
cat(sprintf('SENSITIVITY (re-including day-level misclassifications): HR %.3f (%.3f-%.3f)\n',
            b$estimate, b[["2.5 %"]], b[["97.5 %"]]))
write.csv(data.frame(analysis='Re-including patients with precise deathtime after the study',
                     n=nrow(A2), deaths=sum(A2$event90), hr=b$estimate, lo=b[["2.5 %"]], hi=b[["97.5 %"]]),
          'TableS8b_reinclude.csv', row.names=FALSE)
