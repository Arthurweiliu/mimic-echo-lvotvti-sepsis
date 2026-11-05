## 55_ph_mi.R — proportional-hazards check and episode-specific HRs in the MI cohort
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(mice)); suppressMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(x))

d <- fread(file.path(DATA, 'h1_data.csv'), colClasses=list(character=c('stay_id','subject_id')))
d[, icu_intime := as.POSIXct(icu_intime, tz='UTC')]
for (v in c('age','sofa','sapsii','charlson','lactate','creat','hgb','plt','wbc','bili','inr',
            'hr_mean','mbp_min','temp_mean','spo2_min','bun','norepi','fu_days','event90','vti','vti_diam','lvef'))
  d[, (v) := num(get(v))]
d[, sex := factor(ifelse(sex %in% c('M','F'), sex, NA))]
tf <- function(x) { s <- trimws(as.character(x)); s[is.na(s)] <- ''; ifelse(grepl('^t', tolower(s)), 1, 0) }
d[, `:=`(ventN = tf(vent), rrtN = tf(rrt))]
lk <- fread(file.path(DATA, "linked_cohort.csv"), colClasses=list(character='stay_id'))
lk <- lk[!duplicated(stay_id), .(stay_id, echo_dt_minus_intime_h)]
d <- merge(d, lk, by='stay_id', all.x=TRUE, sort=FALSE)
d[, echo_h := num(echo_dt_minus_intime_h)]
ic <- fread(file.path(DATA, "sepsis3_icustays.csv"), colClasses=list(character=c('subject_id','stay_id')))
ic[, intime := as.POSIXct(intime, tz='UTC')]; setorder(ic, subject_id, intime)
first_ids <- ic[!duplicated(subject_id), stay_id]
d[, start := echo_h/24]
A <- d[!is.na(vti) & vti>0 & !is.na(fu_days) & fu_days>0 & !is.na(echo_h) &
       echo_h>=0 & echo_h<=72 & fu_days>start & stay_id %in% first_ids]
A <- merge(A, fread('lactate_v2.csv', colClasses=list(character=c('stay_id','subject_id'))),
           by='stay_id', all.x=TRUE, suffixes=c('', '.v2'))
A[, vlow := as.numeric(vti < 17)]
aux <- c('sapsii','mbp_min','hr_mean','temp_mean','spo2_min','bun','hgb','plt','wbc','bili',
         'inr','norepi','ventN','rrtN','vti','vti_diam')
D <- unique(A[, c('vlow','age','sex','sofa','charlson','lactate_filled2','creat','event90','fu_days','start', aux), with=FALSE])
sf <- survfit(Surv(start, fu_days, event90) ~ 1, data=D)
D[, nag := approx(sf$time, sf$cumhaz, fu_days, rule=2, ties='ordered')$y]

## Schoenfeld test in the complete-case fit
cc <- D[complete.cases(D[, c('vlow','age','sex','sofa','charlson','lactate_filled2','creat'), with=FALSE])]
f0 <- coxph(Surv(start, fu_days, event90) ~ vlow + age + sex + sofa + charlson + lactate_filled2 + creat, data=cc)
z <- cox.zph(f0)
cat('Schoenfeld (complete case, n =', nrow(cc), '):\n'); print(z$table)

## MI with a 7-day split
M <- as.data.frame(D)
ini <- mice(M, m=1, maxit=0, printFlag=FALSE); meth <- ini$method
for (v in names(meth)) meth[v] <- ''
for (v in c('lactate_filled2','creat', aux)) if (sum(is.na(M[[v]])) > 0) meth[v] <- 'pmm'
set.seed(20261105)
imp <- mice(M, m=20, maxit=10, method=meth, printFlag=FALSE, seed=20261105)

split7 <- function(dk) {
  a <- dk[, c('vlow','age','sex','sofa','charlson','lactate_filled2','creat','start','fu_days','event90')]
  a$t0 <- a$start; a$t1 <- pmin(a$fu_days, 7); a$ev <- ifelse(a$fu_days <= 7, a$event90, 0); a$period <- 'days 0-7'
  b <- a[a$t0 < 7, ]
  c2 <- a[a$fu_days > 7, ]; c2$t0 <- pmax(c2$start, 7); c2$t1 <- c2$fu_days; c2$ev <- c2$event90; c2$period <- 'days 7-90'
  rbind(b, c2)
}
fits <- lapply(seq_len(imp$m), function(k) {
  dk <- split7(complete(imp, k))
  coxph(Surv(t0, t1, ev) ~ vlow*relevel(factor(period), 'days 0-7') + age + sex + sofa + charlson + lactate_filled2 + creat, data=dk)
})
pl <- summary(pool(as.mira(fits)), conf.int=TRUE, exponentiate=TRUE)
print(pl[, c('term','estimate','2.5 %','97.5 %','p.value')])

## day-7 landmark: late (days 7-90) hazard in survivors to day 7
D7 <- D[fu_days > 7, ]; D7[, `:=`(start7=7)]
M7 <- as.data.frame(D7)
ini7 <- mice(M7, m=1, maxit=0, printFlag=FALSE); meth7 <- ini7$method
for (v in names(meth7)) meth7[v] <- ''
for (v in c('lactate_filled2','creat', aux)) if (sum(is.na(M7[[v]])) > 0) meth7[v] <- 'pmm'
set.seed(20261105)
imp7 <- mice(M7, m=20, maxit=10, method=meth7, printFlag=FALSE, seed=20261105)
f7 <- lapply(seq_len(imp7$m), function(k) coxph(Surv(start7, fu_days, event90) ~ vlow + age + sex + sofa +
        charlson + lactate_filled2 + creat, data=complete(imp7, k)))
pl7 <- summary(pool(as.mira(f7)), conf.int=TRUE, exponentiate=TRUE)
b7 <- pl7[grep('^vlow', pl7$term),]
cat(sprintf('\nDay-7 landmark (n=%d, deaths=%d): HR %.3f (%.3f-%.3f)\n', nrow(D7), sum(D7$event90),
            b7$estimate, b7[["2.5 %"]], b7[["97.5 %"]]))
