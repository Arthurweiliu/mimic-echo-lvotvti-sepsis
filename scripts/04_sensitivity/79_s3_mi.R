## 79_s3_mi.R — functional form for every same-study measure under multiple imputation
## For each of VTI, SV, CO, CI, HR, LVEF: impute covariates (m = 20), fit a linear model on the
## untransformed exposure and a 3-df natural spline, pool the linear coefficient (Rubin) and
## report the median likelihood-ratio p for non-linearity and overall association.
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(splines)); suppressMessages(library(mice))
suppressMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(x))
set.seed(20261105)

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
  , .(stay_id, lactate_filled2)], by='stay_id', all.x=TRUE, sort=FALSE)
A[, hrv := ifelse(!is.na(resting_hr) & resting_hr > 25 & resting_hr < 220, resting_hr, NA)]
A[, sv := vti*pi*(vti_diam/2)^2]
A[, co := sv*hrv/1000]
A[, ci := ifelse(!is.na(bsa) & bsa >= 0.9 & bsa <= 3.0 & sv >= 10 & sv <= 250, co/bsa, NA)]

AUX <- c('sapsii','mbp_min','hr_mean','temp_mean','spo2_min','bun','hgb','plt','wbc','bili',
         'inr','norepi','ventN','rrtN','vti_diam')
COV <- c('age','sex','sofa','charlson','lactate_filled2','creat')

analyse <- function(v, unit, label, logx) {
  dat <- A[!is.na(A[[v]]) & A[[v]] > 0]
  ## impute covariates WITHOUT any derived columns of the exposure in the predictor matrix
  X <- as.data.frame(dat[, c('start','fu_days','event90', v, COV, AUX), with=FALSE])
  xraw <- X[[v]]
  xv <- if (logx) log(xraw) else xraw
  B <- ns(xv, df=3); colnames(B) <- paste0('s', 1:3)
  meth <- rep('', ncol(X)); names(meth) <- names(X)
  for (vv in c(COV, AUX)) if (sum(is.na(X[[vv]])) > 0) meth[vv] <- 'pmm'
  sf <- survfit(Surv(start, fu_days, event90) ~ 1, data=X)
  X$nag <- approx(sf$time, sf$cumhaz, X$fu_days, rule=2, ties='ordered')$y
  X$nag[is.na(X$nag)] <- 0
  meth <- c(meth, nag='')
  imp <- mice(X, m=20, maxit=10, method=meth, printFlag=FALSE, seed=20261105)
  get <- function(dk, rhs) coxph(as.formula(paste('Surv(start,fu_days,event90) ~', paste(COV, collapse=' + '), '+', rhs)), data=dk)
  coefs <- ses <- numeric(imp$m); pnl <- pall <- daic <- numeric(imp$m)
  for (k in seq_len(imp$m)) {
    dk <- cbind(complete(imp, k), B)
    dk$lin <- dk[[v]]
    dk$ln  <- log(dk[[v]])
    f0 <- get(dk, '1'); fl <- get(dk, 'lin'); fs <- get(dk, 's1 + s2 + s3')
    f1 <- if (logx) get(dk, 'ln') else fl
    coefs[k] <- coef(fl)['lin']; ses[k] <- sqrt(vcov(fl)['lin','lin'])
    pnl[k]  <- pchisq(2*(fs$loglik[2] - f1$loglik[2]), 2, lower.tail=FALSE)
    pall[k] <- pchisq(2*(fs$loglik[2] - f0$loglik[2]), 3, lower.tail=FALSE)
    daic[k] <- AIC(f1) - AIC(fs)
  }
  m <- imp$m; b <- mean(coefs); T <- mean(ses^2) + (1 + 1/m)*var(coefs); se <- sqrt(T)
  data.frame(measure=label, unit=sprintf('per %g', unit), n=nrow(X), deaths=sum(X$event90),
             hr=exp(b*unit), lo=exp((b-1.96*se)*unit), hi=exp((b+1.96*se)*unit),
             p_lin=2*pnorm(-abs(b/se)), p_nonlin=median(pnl), p_overall=median(pall),
             dAIC=median(daic), pct_nl_sig=100*mean(pnl < 0.05))
}

res <- rbind(
  analyse('vti', 5, 'LVOT-VTI', TRUE),
  analyse('sv', 10, 'Stroke volume', TRUE),
  analyse('co', 1, 'Cardiac output', TRUE),
  analyse('ci', 1, 'Cardiac index', TRUE),
  analyse('hrv', 10, 'Heart rate', FALSE),
  analyse('lvef', 10, 'LVEF', FALSE))
print(res, digits=3)
write.csv(res, 'TableS3_mi.csv', row.names=FALSE)

fmt <- function(x, d) formatC(x, format='f', digits=d)
md <- c('### Table S3. Linear versus non-linear functional form for every same-study measure (multiple imputation)', '',
        '| Measure | Contrast | n | Deaths | Linear HR (95% CI) | p linear | Non-linearity p | Overall p | delta AIC (spline - linear) |',
        '|---|---|---:|---:|---|---:|---:|---:|---:|')
for (i in seq_len(nrow(res))) with(res[i, ], md <<- c(md, sprintf('| %s | %s | %s | %s | %s (%s-%s) | %s | %s | %s | %+.1f |',
        measure, unit, format(n, big.mark=','), format(deaths, big.mark=','),
        fmt(hr, 2), fmt(lo, 2), fmt(hi, 2),
        if (p_lin < 0.001) '<0.001' else fmt(p_lin, 3),
        formatC(p_nonlin, format='g', digits=2), formatC(p_overall, format='g', digits=2), dAIC)))
md <- c(md, '', paste0('*Covariates imputed (m = 20) within each measure-specific cohort; the markers themselves are not imputed. The linear coefficient is fitted on the untransformed exposure (the quantity under test) and pooled by Rubin\'s rules; the non-linearity test compares a 3-df natural spline with a linear term on the same scale (log for the flow measures, untransformed for LVEF and heart rate) and is reported as the median p across the 20 imputations. A positive delta AIC means the spline fits better.*'))
writeLines(md, 'TableS3_mi.md')
cat('\nwritten TableS3_mi.csv/.md\n')
