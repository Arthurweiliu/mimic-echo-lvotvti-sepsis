## 67_mi_headtohead.R — head-to-head models in the multiply imputed cohort
## H1: all patients with VTI + LVEF from the same index TTE
## H2: those with a computable echo-derived CI
## Covariates are imputed (m = 20); the markers are never imputed.
## Predicted 90-day risks are averaged over imputations; discrimination and calibration
## are then computed on the averaged risk.
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(mice)); suppressMessages(library(splines))
suppressMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(x))
set.seed(20261105)

## ---------- cohort ----------
d <- fread(file.path(DATA, "h1_data.csv"), colClasses=list(character=c('stay_id','subject_id')))
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
A <- merge(A, fread(file.path(DATA, "analysis_data.csv"),
                    colClasses=list(character=c('stay_id','measurement_id')))[
  !duplicated(stay_id), .(stay_id, measurement_id=as.character(measurement_id), bsa=num(body_surface_area))],
  by='stay_id', all.x=TRUE, sort=FALSE)
ex <- read.csv(file.path(DATA, "echo_extra.csv"), stringsAsFactors=FALSE,
               colClasses=c(measurement_id='character'))
hr <- data.table(measurement_id=ex$measurement_id[ex$measurement=='resting_hr'],
                 resting_hr=num(ex$result[ex$measurement=='resting_hr']))[!duplicated(measurement_id)]
A <- merge(A, hr, by='measurement_id', all.x=TRUE, sort=FALSE)
A <- merge(A, fread('lactate_v2.csv', colClasses=list(character=c('stay_id','subject_id')))[
  , .(stay_id, lactate_filled2, lactate_nearest)], by='stay_id', all.x=TRUE, sort=FALSE)
A[, vlow := as.numeric(vti < 17)]
A[, hrv := ifelse(!is.na(resting_hr) & resting_hr > 25 & resting_hr < 220, resting_hr, NA)]
A[, sv := vti*pi*(vti_diam/2)^2]
A[, co := sv*hrv/1000]
A[, ok_ci := vti >= 5 & vti <= 45 & !is.na(vti_diam) & vti_diam >= 1.0 & vti_diam <= 3.2 &
    !is.na(hrv) & !is.na(bsa) & bsa >= 0.9 & bsa <= 3.0 & sv >= 10 & sv <= 250]
A[, ci := ifelse(ok_ci, co/bsa, NA)]
cat('cohort', nrow(A), '| derived-CI plausibility exclusions:', sum(!A$ok_ci), '\n')

AUX <- c('sapsii','mbp_min','hr_mean','temp_mean','spo2_min','bun','hgb','plt','wbc','bili',
         'inr','norepi','ventN','rrtN','vti_diam')
COV <- c('age','sex','sofa','charlson','lactate_filled2','creat')

mk_frame <- function(dat, with_eci) {
  X <- as.data.frame(dat[, c('start','fu_days','event90','vti','lvef', if (with_eci) 'ci', 'vlow', COV, AUX), with=FALSE])
  B <- ns(log(X$vti), df=3); colnames(B) <- paste0('v_s', 1:3); X <- cbind(X, B)
  B <- ns(X$lvef, df=3); colnames(B) <- paste0('l_s', 1:3); X <- cbind(X, B)
  if (with_eci) { B <- ns(log(X$ci), df=3); colnames(B) <- paste0('e_s', 1:3); X <- cbind(X, B) }
  X
}

auc <- function(y, p) { r <- rank(p); n1 <- sum(y==1); n0 <- sum(y==0); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }
## survival::concordance() returns 1 - C for a risk score under left truncation in this build;
## the formula form with reverse=TRUE reproduces a manual Harrell C (verified: 0.6991 both)
cindex <- function(X, p) as.numeric(concordance(Surv(X$start, X$fu_days, X$event90) ~ p, reverse=TRUE)$concordance)
brier <- function(X, p) mean((X$event90 - p)^2)
calib <- function(X, p) { m <- glm(X$event90 ~ qlogis(p), family=binomial)
  c(slope=unname(coef(m)[2]), intercept=unname(coef(m)[1])) }
hl <- function(X, p, g=10) {
  br <- unique(quantile(p, seq(0, 1, 1/g))); q <- cut(p, br, include.lowest=TRUE)
  obs <- tapply(X$event90, q, sum); ex <- tapply(p, q, sum)
  stat <- sum((obs - ex)^2 / pmax(ex, 1e-9)); df <- length(obs) - 2
  c(chisq=stat, df=df, p=pchisq(stat, df, lower.tail=FALSE))
}

run_set <- function(dat, models, label) {
  X <- mk_frame(dat, with_eci=any(grepl('e_s', models)))
  cat(sprintf('\n===== %s: n = %d, deaths = %d =====\n', label, nrow(X), sum(X$event90)))
  meth <- rep('', ncol(X)); names(meth) <- names(X)
  for (v in c(COV, AUX)) if (v %in% names(X) && sum(is.na(X[[v]])) > 0) meth[v] <- 'pmm'
  sf <- survfit(Surv(start, fu_days, event90) ~ 1, data=X)
  X$nag <- approx(sf$time, sf$cumhaz, X$fu_days, rule=2, ties='ordered')$y
  X$nag[is.na(X$nag)] <- 0
  meth <- c(meth, nag = '')
  imp <- mice(X, m=20, maxit=10, method=meth, printFlag=FALSE, seed=20261105)
  comp <- lapply(seq_len(imp$m), function(k) complete(imp, k))
  risks <- array(NA_real_, dim=c(nrow(X), imp$m, length(models)),
                 dimnames=list(NULL, NULL, names(models)))
  for (k in seq_len(imp$m)) {
    dk <- comp[[k]]
    for (nm in names(models)) {
      f <- coxph(as.formula(models[[nm]]), data=dk)
      bh <- basehaz(f, centered=FALSE)
      H90 <- approx(bh$time, bh$hazard, 90, rule=2)$y
      Hs <- approx(bh$time, bh$hazard, dk$start, rule=2)$y
      lp <- predict(f, newdata=dk, type='lp', reference='zero')
      risks[, k, nm] <- pmin(pmax(1 - exp(-(H90 - Hs) * exp(lp)), 1e-6), 1 - 1e-6)
    }
  }
  avg <- apply(risks, c(1, 3), mean)
  y <- X$event90
  base <- 'clinical'
  tab <- do.call(rbind, lapply(names(models), function(nm) {
    p <- avg[, nm]; cb <- calib(X, p); h <- hl(X, p)
    data.frame(set=label, model=nm, n=nrow(X), deaths=sum(y), AUC=auc(y, p),
               Cindex=cindex(X, p), Brier=brier(X, p), slope=cb['slope'],
               intercept=cb['intercept'], HL_p=h['p'], dAUC=auc(y, p) - auc(y, avg[, base]))
  }))
  print(tab, digits=4)
  list(X=X, avg=avg, comp=comp, tab=tab, imp=imp)
}

H1 <- A[!is.na(lvef) & lvef > 0]
H2 <- H1[!is.na(ci)]
cat('H1 n', nrow(H1), 'deaths', sum(H1$event90), '| H2 n', nrow(H2), 'deaths', sum(H2$event90), '\n')

covs <- 'age + sex + sofa + charlson + lactate_filled2 + creat'
M1 <- list(clinical=paste('Surv(start,fu_days,event90) ~', covs),
           VTI=paste('Surv(start,fu_days,event90) ~', covs, '+ v_s1 + v_s2 + v_s3'),
           LVEF=paste('Surv(start,fu_days,event90) ~', covs, '+ l_s1 + l_s2 + l_s3'),
           VTI_LVEF=paste('Surv(start,fu_days,event90) ~', covs, '+ v_s1 + v_s2 + v_s3 + l_s1 + l_s2 + l_s3'))
M2 <- c(M1, list(
  eCI=paste('Surv(start,fu_days,event90) ~', covs, '+ e_s1 + e_s2 + e_s3'),
  VTI_eCI=paste('Surv(start,fu_days,event90) ~', covs, '+ v_s1 + v_s2 + v_s3 + e_s1 + e_s2 + e_s3'),
  LVEF_eCI=paste('Surv(start,fu_days,event90) ~', covs, '+ l_s1 + l_s2 + l_s3 + e_s1 + e_s2 + e_s3'),
  all=paste('Surv(start,fu_days,event90) ~', covs,
            '+ v_s1 + v_s2 + v_s3 + l_s1 + l_s2 + l_s3 + e_s1 + e_s2 + e_s3')))

r1 <- run_set(H1, M1, 'H1')
r2 <- run_set(H2, M2, 'H2')
saveRDS(list(H1=r1$tab, H2=r2$tab), 'Table3_mi.rds')
write.csv(rbind(r1$tab, r2$tab), 'Table3_mi.csv', row.names=FALSE)
saveRDS(list(X=r1$X, comp=r1$comp, avg=r1$avg), 'mi_comp_H1.rds')
saveRDS(list(X=r2$X, comp=r2$comp, avg=r2$avg), 'mi_comp_H2.rds')

## ---------- per-SD effects, mutually adjusted H2 linear model (MI) ----------
X2 <- r2$X
sdv <- c(VTI=sd(log(X2$vti)), LVEF=sd(X2$lvef), eCI=sd(log(X2$ci)))  # log scale for log-modelled terms
cat('\nSDs on the modelled scale:', round(sdv, 3), '\n')
fitps <- lapply(r2$comp, function(dk)
  coxph(Surv(start,fu_days,event90) ~ age + sex + sofa + charlson + lactate_filled2 + creat +
          log(vti) + lvef + log(ci), data=dk))
pool_ps <- function(term, mult) {
  bs <- sapply(fitps, function(f) coef(f)[term]); ses <- sapply(fitps, function(f) sqrt(vcov(f)[term, term]))
  m <- mean(bs); T <- mean(ses^2) + (1 + 1/length(bs))*var(bs); se <- sqrt(T)
  c(est=exp(m*mult), lo=exp((m - 1.96*se)*mult), hi=exp((m + 1.96*se)*mult), p=2*pnorm(-abs(m/se)))
}
ps <- rbind(VTI=pool_ps('log(vti)', sdv['VTI']), LVEF=pool_ps('lvef', sdv['LVEF']),
            eCI=pool_ps('log(ci)', sdv['eCI']))
print(round(ps, 3))
write.csv(ps, 'Table4_mi_perSD.csv')
cat('\nwritten Table3_mi.csv, Table4_mi_perSD.csv, mi_comp_H1/H2.rds\n')
