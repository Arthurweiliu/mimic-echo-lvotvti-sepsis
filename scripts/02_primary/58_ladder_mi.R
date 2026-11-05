## 58_ladder_mi.R — nested adjustment ladder in the MI cohort (n = 4,816)
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(mice)); suppressMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(x))
d <- fread(file.path(DATA, 'h1_data.csv'), colClasses=list(character=c('stay_id','subject_id')))
d[, icu_intime := as.POSIXct(icu_intime, tz='UTC')]
for (v in c('age','sofa','sapsii','charlson','creat','hgb','plt','wbc','bili','inr','hr_mean',
            'mbp_min','temp_mean','spo2_min','bun','norepi','fu_days','event90','vti','vti_diam','lvef'))
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
A <- d[!is.na(vti) & vti>0 & !is.na(fu_days) & fu_days>0 & !is.na(echo_h) & echo_h>=0 &
       echo_h<=72 & fu_days>start & stay_id %in% first_ids]
A <- merge(A, fread('lactate_v2.csv', colClasses=list(character=c('stay_id','subject_id'))),
           by='stay_id', all.x=TRUE, suffixes=c('', '.v2'))
A[, vlow := as.numeric(vti < 17)]
aux <- c('sapsii','mbp_min','hr_mean','temp_mean','spo2_min','bun','hgb','plt','wbc','bili',
         'inr','norepi','ventN','rrtN','vti_diam')
D <- as.data.frame(A[, c('vlow','age','sex','sofa','charlson','lactate_filled2','creat','event90',
                         'fu_days','start', aux), with=FALSE])
sf <- survfit(Surv(start, fu_days, event90) ~ 1, data=D)
D$nag <- approx(sf$time, sf$cumhaz, D$fu_days, rule=2, ties='ordered')$y
ini <- mice(D, m=1, maxit=0, printFlag=FALSE); meth <- ini$method
for (v in names(meth)) meth[v] <- ''
for (v in c('lactate_filled2','creat', aux)) if (sum(is.na(D[[v]])) > 0) meth[v] <- 'pmm'
set.seed(20261105)
imp <- mice(D, m=20, maxit=10, method=meth, printFlag=FALSE, seed=20261105)

show <- function(form, label) {
  fits <- lapply(seq_len(imp$m), function(k) coxph(form, data=complete(imp, k)))
  s <- summary(pool(as.mira(fits)), conf.int=TRUE, exponentiate=TRUE)
  b <- s[grep('^vlow', s$term),]
  cat(sprintf('%-42s HR %.3f (%.3f-%.3f)\n', label, b$estimate, b[["2.5 %"]], b[["97.5 %"]]))
}
show(Surv(start, fu_days, event90) ~ vlow, 'crude')
show(Surv(start, fu_days, event90) ~ vlow + age + sex + charlson, '+ age, sex, Charlson')
show(Surv(start, fu_days, event90) ~ vlow + age + sex + charlson + sofa + lactate_filled2 + creat + norepi + ventN, '+ SOFA, lactate, creatinine, vasopressor, ventilation')
