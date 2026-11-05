## 52_mi_primary_v2.R — PRIMARY analysis: MI (m=20) on the full eligible cohort
## Lactate: blood-gas first-day value; supplemented with chemistry lactate where missing.
## Produces the missing-data cascade table used in the manuscript.
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(mice)); suppressMessages(library(data.table))

num <- function(x) suppressWarnings(as.numeric(x))
d <- fread(file.path(DATA, 'h1_data.csv'), colClasses=list(character=c('stay_id','subject_id')))
d[, icu_intime := as.POSIXct(icu_intime, tz='UTC')]
for (v in c('age','sofa','sapsii','charlson','lactate','creat','hgb','plt','wbc','bili','inr',
            'hr_mean','mbp_min','temp_mean','spo2_min','bun','norepi','fu_days','event90',
            'vti','vti_diam','lvef','troponin'))
  d[, (v) := num(get(v))]
d[, sex := factor(ifelse(sex %in% c('M','F'), sex, NA))]
tf <- function(x) ifelse(grepl('^t', tolower(trimws(ifelse(is.na(x), '', as.character(x))))), 1,
                  ifelse(trimws(ifelse(is.na(x), '', as.character(x))) == 'f', 0, NA))
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

cat('=== Full eligible cohort ===\n')
cat('n =', nrow(A), '| deaths =', sum(A$event90), '\n')
cat('lactate source:\n'); print(A[, .N, by=lactate_src])
cat('missing lactate (bg+chemistry):', sum(is.na(A$lactate_filled2)),
    sprintf('(%.1f%%)', 100*mean(is.na(A$lactate_filled2))),
    '| after nearest rule:', sum(is.na(A$lactate_nearest)),
    sprintf('(%.1f%%)', 100*mean(is.na(A$lactate_nearest))), '\n')

## ---- imputation frame -------------------------------------------------
aux <- c('sapsii','mbp_min','hr_mean','temp_mean','spo2_min','bun','hgb','plt','wbc','bili',
         'inr','norepi','ventN','rrtN','vti','vti_diam')
D <- A[, c('vlow','age','sex','sofa','charlson','lactate','lactate_filled2','lactate_nearest',
           'creat','event90','fu_days','start','vti', aux), with=FALSE]
sf <- survfit(Surv(start, fu_days, event90) ~ 1, data=D)
D[, nag := approx(sf$time, sf$cumhaz, fu_days, rule=2, ties='ordered')$y]

mk_meth <- function(Mf, targets) {
  ini <- mice(Mf, m=1, maxit=0, printFlag=FALSE); mt <- ini$method
  for (v in names(mt)) mt[v] <- ''
  for (v in c(targets, aux)) if (v %in% names(mt) && sum(is.na(Mf[[v]])) > 0) mt[v] <- 'pmm'
  mt
}
run_mi <- function(dat, lacvar, label, m=20) {
  keep <- c('vlow','age','sex','sofa','charlson', lacvar, 'creat','event90','nag','start','fu_days', aux)
  M <- as.data.frame(dat[, keep, with=FALSE])
  set.seed(20261105)
  im <- mice(M, m=m, maxit=10, method=mk_meth(M, lacvar), printFlag=FALSE, seed=20261105)
  if (!is.null(im$loggedEvents)) print(im$loggedEvents)
  form <- as.formula(paste('Surv(start, fu_days, event90) ~ vlow + age + sex + sofa + charlson +',
                           lacvar, '+ creat'))
  ft <- pool(as.mira(lapply(seq_len(im$m), function(k) coxph(form, data=complete(im, k)))))
  pl <- ft; s <- summary(pl, conf.int=TRUE, exponentiate=TRUE)
  b <- s[grep('^vlow', s$term),]
  fmi <- pl$pool$fmi[match(b$term, pl$pool$term)]
  cat(sprintf('%-46s n=%d deaths=%d | HR %.3f (%.3f-%.3f) p=%.2g FMI=%.3f\n', label,
              nrow(M), sum(M$event90), b$estimate, b[["2.5 %"]], b[["97.5 %"]], b$p.value, fmi))
  invisible(list(s=s, imp=im))
}
cat('\n--- MI models ---\n')
r_bg  <- run_mi(D, 'lactate',          'MI, blood-gas lactate only')
r_fil <- run_mi(D, 'lactate_filled2',  'MI, blood gas + chemistry in window (PRIMARY)')
r_nea <- run_mi(D, 'lactate_nearest',  'MI, blood gas + chemistry + nearest (sens.)')

## pooled estimates for all covariates of the primary model
s <- r_fil$s
cat('\nPrimary model, pooled coefficients:\n')
print(s[, c('term','estimate','2.5 %','97.5 %','p.value')])

## per +5 cm and <17 vs >20 in the primary imputation
M <- D[, c('vlow','age','sex','sofa','charlson','lactate_filled2','creat','event90','nag','start','fu_days', aux), with=FALSE]
set.seed(20261105)
imp <- r_fil$imp
f5 <- pool(as.mira(lapply(seq_len(imp$m), function(k) coxph(Surv(start, fu_days, event90) ~ I(vti/5) + age + sex + sofa + charlson + lactate_filled2 + creat, data=complete(imp,k)))))
s5 <- summary(f5, conf.int=TRUE, exponentiate=TRUE); b5 <- s5[grep('vti/5', s5$term),]
cat(sprintf('\nper +5 cm: %.3f (%.3f-%.3f)\n', b5$estimate, b5[["2.5 %"]], b5[["97.5 %"]]))

Dsub <- D[vti<17 | vti>20,]; Dsub[, vlow20 := as.numeric(vti<17)]
set.seed(20261105)
Ms <- as.data.frame(Dsub[, c('vlow20','age','sex','sofa','charlson','lactate_filled2','creat','event90','nag','start','fu_days', aux), with=FALSE])
ims <- mice(Ms, m=20, maxit=10, method=mk_meth(Ms, 'lactate_filled2'), printFlag=FALSE, seed=20261105)
fs <- pool(as.mira(lapply(seq_len(ims$m), function(k) coxph(Surv(start, fu_days, event90) ~ vlow20 + age + sex + sofa + charlson + lactate_filled2 + creat, data=complete(ims,k)))))
ss <- summary(fs, conf.int=TRUE, exponentiate=TRUE); bs <- ss[grep('vlow20', ss$term),]
cat(sprintf('<17 vs >20 (n=%d): %.3f (%.3f-%.3f)\n', nrow(Dsub), bs$estimate, bs[["2.5 %"]], bs[["97.5 %"]]))

## 72-h landmark
Dl <- D[fu_days > 3,]; Dl[, `:=`(startL=3, vlowL=vlow)]
set.seed(20261105)
Ml <- as.data.frame(Dl[, c('vlowL','age','sex','sofa','charlson','lactate_filled2','creat','event90','nag','startL','fu_days', aux), with=FALSE])
ml <- mice(Ml, m=20, maxit=10, method=mk_meth(Ml, 'lactate_filled2'), printFlag=FALSE, seed=20261105)
fl <- pool(as.mira(lapply(seq_len(ml$m), function(k) coxph(Surv(startL, fu_days, event90) ~ vlowL + age + sex + sofa + charlson + lactate_filled2 + creat, data=complete(ml,k)))))
sl <- summary(fl, conf.int=TRUE, exponentiate=TRUE); bl <- sl[grep('vlowL', sl$term),]
cat(sprintf('72-h landmark (n=%d, deaths=%d): %.3f (%.3f-%.3f)\n', nrow(Dl), sum(Dl$event90),
            bl$estimate, bl[["2.5 %"]], bl[["97.5 %"]]))

## ---- complete-case rows ----------------------------------------------
cc_row <- function(dat, lacvar, label) {
  ok <- complete.cases(dat[, c('vlow','age','sex','sofa','charlson', lacvar, 'creat'), with=FALSE])
  x <- dat[ok,]
  f <- coxph(as.formula(paste('Surv(start, fu_days, event90) ~ vlow + age + sex + sofa + charlson +',
                              lacvar, '+ creat')), data=x)
  b <- summary(f)$coefficients['vlow',]
  cat(sprintf('%-46s n=%d deaths=%d | HR %.3f (%.3f-%.3f)\n', label, nrow(x), sum(x$event90),
              exp(b[1]), exp(b[1]-1.96*b[3]), exp(b[1]+1.96*b[3])))
  invisible(x)
}
cat('\n--- Complete-case models ---\n')
cc_row(D, 'lactate',          'CC, blood-gas lactate only')
cc_row(D, 'lactate_filled2',  'CC, blood gas + chemistry in window')
cc_row(D, 'lactate_nearest',  'CC, blood gas + chemistry + nearest')

## ---- H1 / H2 with MI -------------------------------------------------
cat('\n--- H1 / H2 (selection by same-TTE availability, covariates MI) ---\n')
H1 <- A[!is.na(lvef) & lvef>0,]
cat('H1 n =', nrow(H1), 'deaths =', sum(H1$event90), '\n')
H1d <- as.data.frame(H1[, c('vlow','age','sex','sofa','charlson','lactate_filled2','creat','event90','fu_days','start','vti', aux), with=FALSE])
sf1 <- survfit(Surv(start, fu_days, event90) ~ 1, data=H1d)
H1d$nag <- approx(sf1$time, sf1$cumhaz, H1d$fu_days, rule=2, ties='ordered')$y
set.seed(20261105)
im1 <- mice(H1d, m=20, maxit=10, method=mk_meth(H1d, 'lactate_filled2'), printFlag=FALSE, seed=20261105)
f1 <- pool(as.mira(lapply(seq_len(im1$m), function(k) coxph(Surv(start, fu_days, event90) ~ vlow + age + sex + sofa + charlson + lactate_filled2 + creat, data=complete(im1,k)))))
s1 <- summary(f1, conf.int=TRUE, exponentiate=TRUE); b1 <- s1[grep('vlow', s1$term),]
cat(sprintf('H1 MI: HR %.3f (%.3f-%.3f)\n', b1$estimate, b1[["2.5 %"]], b1[["97.5 %"]]))

saveRDS(list(D=D, A=A, s=s), 'mi_primary_v2.rds')
cat('\nsaved mi_primary_v2.rds\n')
