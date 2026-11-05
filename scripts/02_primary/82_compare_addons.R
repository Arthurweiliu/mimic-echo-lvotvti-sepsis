## 82_compare_addons.R — two checks suggested by the published LVEF paper
## (1) 28-day outcome for the primary contrast (the field's usual horizon)
## (2) pooled non-linearity LRT across imputations (Meng-Rubin D3) instead of a median p
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(splines)); suppressMessages(library(mice))
suppressMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(x)); set.seed(20261105)

d <- fread(file.path(DATA, "h1_data.csv"), colClasses=list(character=c('stay_id','subject_id')))
for (v in c('age','sofa','charlson','creat','fu_days','event90','death28','vti','lvef','vti_diam'))
  d[, (v) := num(get(v))]
d[, sex := factor(ifelse(sex %in% c('M','F'), sex, NA))]
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
A[, vlow := as.numeric(vti < 17)]
cat('cohort', nrow(A), '| 28-day deaths', sum(A$death28 == 1, na.rm=TRUE),
    '| 90-day deaths', sum(A$event90), '\n')

AUX <- c('sapsii','mbp_min','hr_mean','temp_mean','spo2_min','bun','hgb','plt','wbc','bili',
         'inr','norepi','ventN','rrtN','vti_diam')
tf <- function(x) { s <- trimws(as.character(x)); s[is.na(s)] <- ''; ifelse(grepl('^t', tolower(s)), 1, 0) }
A[, `:=`(ventN = tf(vent), rrtN = tf(rrt))]
COV <- c('age','sex','sofa','charlson','lactate_filled2','creat')

## ---------- (1) 28-day ----------
X <- as.data.frame(A[, c('start','fu_days','event90','death28','vlow','vti', COV, AUX), with=FALSE])
X$t28 <- pmin(X$fu_days, 28)
X$e28 <- as.numeric(X$death28 == 1 & X$fu_days <= 28)
meth <- rep('', ncol(X)); names(meth) <- names(X)
for (v in c(COV, AUX)) if (sum(is.na(X[[v]])) > 0) meth[v] <- 'pmm'
sf <- survfit(Surv(start, fu_days, event90) ~ 1, data=X)
X$nag <- approx(sf$time, sf$cumhaz, X$fu_days, rule=2, ties='ordered')$y
X$nag[is.na(X$nag)] <- 0
meth <- c(meth, nag='')
imp <- mice(X, m=20, maxit=10, method=meth, printFlag=FALSE, seed=20261105)
f28 <- lapply(seq_len(imp$m), function(k)
  coxph(Surv(start, t28, e28) ~ vlow + age + sex + sofa + charlson + lactate_filled2 + creat,
        data=complete(imp, k)))
p28 <- summary(pool(as.mira(f28)), conf.int=TRUE, exponentiate=TRUE)
b <- p28[grep('^vlow', p28$term),]
cat(sprintf('28-day: VTI<17 vs >=17  HR %.3f (%.3f-%.3f)\n', b$estimate, b[["2.5 %"]], b[["97.5 %"]]))

## ---------- (2) pooled D3 for non-linearity ----------
D3test <- function(v, logx, label) {
  dd <- A[!is.na(A[[v]]) & A[[v]] > 0]
  Xv <- as.data.frame(dd[, c('start','fu_days','event90', v, COV, AUX), with=FALSE])
  xv <- if (logx) log(Xv[[v]]) else Xv[[v]]
  B <- ns(xv, df=3); colnames(B) <- paste0('s', 1:3)
  methv <- rep('', ncol(Xv)); names(methv) <- names(Xv)
  for (vv in c(COV, AUX)) if (sum(is.na(Xv[[vv]])) > 0) methv[vv] <- 'pmm'
  sfv <- survfit(Surv(start, fu_days, event90) ~ 1, data=Xv)
  Xv$nag <- approx(sfv$time, sfv$cumhaz, Xv$fu_days, rule=2, ties='ordered')$y
  Xv$nag[is.na(Xv$nag)] <- 0
  methv <- c(methv, nag='')
  imv <- mice(Xv, m=20, maxit=10, method=methv, printFlag=FALSE, seed=20261105)
  comp <- lapply(seq_len(imv$m), function(k) { dk <- cbind(complete(imv,k), B); dk$lx <- log(dk[[v]]); dk$lin <- dk[[v]]; dk })
  base <- paste(COV, collapse=' + ')
  fl <- lapply(comp, function(dk) coxph(as.formula(paste('Surv(start,fu_days,event90) ~', base, '+', if (logx) 'lx' else 'lin')), data=dk))
  fs <- lapply(comp, function(dk) coxph(as.formula(paste('Surv(start,fu_days,event90) ~', base, '+ s1 + s2 + s3')), data=dk))
  d3 <- tryCatch(D3(as.mira(fl), as.mira(fs)), error=function(e) NULL)
  if (is.null(d3)) { cat(label, ': D3 failed\n'); return(invisible()) }
  cat(sprintf('%-6s pooled D3 = %.2f, df = %d, p = %.3g\n', label, d3$D3, d3$df, d3$p))
}
D3test('vti', TRUE, 'VTI')
D3test('lvef', FALSE, 'LVEF')
