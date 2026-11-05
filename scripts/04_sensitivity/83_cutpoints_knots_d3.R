## 83_cutpoints_knots_d3.R — one imputation set, then:
##   (a) alternative VTI cut-points
##   (b) 3-df vs 4-knot (5/35/65/95) spline for VTI
##   (c) pooled non-linearity test (Meng-Rubin D3) for VTI, LVEF, eCI
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(splines)); suppressMessages(library(mice))
suppressMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(x)); set.seed(20261105)

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
A <- merge(A, fread('lactate_v2.csv', colClasses=list(character=c('stay_id','subject_id')))[
  , .(stay_id, lactate_filled2)], by='stay_id', all.x=TRUE, sort=FALSE)
AUX <- c('sapsii','mbp_min','hr_mean','temp_mean','spo2_min','bun','hgb','plt','wbc','bili',
         'inr','norepi','ventN','rrtN','vti_diam')
COV <- c('age','sex','sofa','charlson','lactate_filled2','creat')
cat('cohort', nrow(A), 'deaths', sum(A$event90), '\n')

X <- as.data.frame(A[, c('start','fu_days','event90','vti','lvef', COV, AUX), with=FALSE])
sf <- survfit(Surv(start, fu_days, event90) ~ 1, data=X)
X$nag <- approx(sf$time, sf$cumhaz, X$fu_days, rule=2, ties='ordered')$y
X$nag[is.na(X$nag)] <- 0
## derived columns are carried but never used as predictors
X$lx <- log(X$vti)
B3 <- ns(log(X$vti), df=3); colnames(B3) <- paste0('s', 1:3); X <- cbind(X, B3)
B4 <- ns(log(X$vti), knots=quantile(log(X$vti), c(.05,.35,.65,.95)), intercept=FALSE)
colnames(B4) <- paste0('k', 1:ncol(B4)); X <- cbind(X, B4)
meth <- rep('', ncol(X)); names(meth) <- names(X)
for (v in c(COV, AUX)) if (sum(is.na(X[[v]])) > 0) meth[v] <- 'pmm'
pred <- matrix(1, ncol(X), ncol(X), dimnames=list(names(X), names(X)))
diag(pred) <- 0
pred[, c('lx', paste0('s', 1:3), paste0('k', 1:ncol(B4)))] <- 0   # keep them out of the imputation model
imp <- mice(X, m=20, maxit=10, method=meth, predictorMatrix=pred, printFlag=FALSE, seed=20261105)
base <- paste(COV, collapse=' + ')
fitall <- function(rhs) with(imp, coxph(as.formula(paste('Surv(start,fu_days,event90) ~', base, '+', rhs))))
hr <- function(m, term) { s <- summary(pool(m), conf.int=TRUE, exponentiate=TRUE); b <- s[grep(term, s$term),]
  sprintf('%.3f (%.3f-%.3f)', b$estimate, b[["2.5 %"]], b[["97.5 %"]]) }

## ---------- (a) cut-points ----------
cuts <- list(c(17, 'strict'), c(17, 'inclusive'), c(16.5, 'strict'), c(18, 'strict'), c(16, 'strict'))
cat('\n--- VTI cut-points (same imputations) ---\n')
rows <- list()
for (cc in cuts) {
  th <- as.numeric(cc[1]); inc <- cc[2] == 'inclusive'
  vname <- paste0('c_', gsub('\\.', '', as.character(th)), if (inc) 'inc' else '')
  X2 <- X; X2[[vname]] <- as.numeric(if (inc) X$vti <= th else X$vti < th)
  ## refit with the new column by completing and appending
  f <- lapply(seq_len(imp$m), function(k) { dk <- cbind(complete(imp, k), setNames(data.frame(X2[[vname]]), vname))
    coxph(as.formula(paste('Surv(start,fu_days,event90) ~', vname, '+', base)), data=dk) })
  m <- as.mira(f)
  nlow <- sum(X2[[vname]] == 1)
  cat(sprintf('  %-22s n<%g = %4d | HR %s\n', paste0(if (inc) '<=' else '<', th, ' vs rest'), th, nlow, hr(m, vname)))
  rows[[length(rows)+1]] <- data.frame(definition = paste0(if (inc) '<=' else '<', th),
                                       n_low = nlow, n_high = nrow(X2) - nlow, hr = hr(m, vname))
}
write.csv(do.call(rbind, rows), 'TableS6_cutpoints.csv', row.names=FALSE)

## ---------- (b) 3-df vs 4-knot spline ----------
cat('\n--- spline specification ---\n')
m3 <- fitall('s1 + s2 + s3')
m4 <- fitall('k1 + k2 + k3 + k4')
lik <- function(m) sapply(m$analyses, function(f) f$loglik[2])
cat(sprintf('  total log-lik (3 df): %.1f | (4 knots): %.1f | dAIC %.1f\n',
            sum(lik(m3)) - 20*0, sum(lik(m4)), 2*(sum(lik(m4)) - sum(lik(m3))) - 2*(4-3)))
d3 <- tryCatch(D3(m3, m4), error=function(e) e)
if (inherits(d3, 'error')) cat('  D3(3 vs 4 df) failed:', conditionMessage(d3), '\n') else
  cat(sprintf('  D3 = %.2f, df = %d, p = %.3g\n', d3$D3, d3$df, d3$p))

## ---------- (c) pooled non-linearity tests ----------
cat('\n--- pooled non-linearity (D3) ---\n')
m_lin <- fitall('lx'); m_spl <- fitall('s1 + s2 + s3')
d <- tryCatch(D3(m_lin, m_spl), error=function(e) e)
if (inherits(d, 'error')) cat('  VTI: D3 failed:', conditionMessage(d), '\n') else
  cat(sprintf('  VTI : D3 = %.2f, df = %d, p = %.3g\n', d$D3, d$df, d$p))
