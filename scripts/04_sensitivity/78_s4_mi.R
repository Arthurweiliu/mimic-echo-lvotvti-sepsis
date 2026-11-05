## 78_s4_mi.R — residual sensitivity analyses under multiple imputation (Table S4, MI)
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(mice)); suppressMessages(library(data.table))
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
d <- merge(d, fread('lactate_v2.csv', colClasses=list(character=c('stay_id','subject_id')))[
  , .(stay_id, lactate_filled2)], by='stay_id', all.x=TRUE, sort=FALSE)

inwin <- !is.na(d$vti) & d$vti > 0 & !is.na(d$echo_h) & d$echo_h >= 0 & d$echo_h <= 72
alive <- inwin & !is.na(d$fu_days) & d$fu_days > d$start
d[, vlow := as.numeric(vti < 17)]
A <- d[alive & d$stay_id %in% first_ids]

AUTO <- c('sapsii','mbp_min','hr_mean','temp_mean','spo2_min','bun','hgb','plt','wbc','bili',
          'inr','norepi','ventN','rrtN','vti_diam')
COV <- c('age','sex','sofa','charlson','lactate_filled2','creat')

run <- function(dat, label) {
  X <- as.data.frame(dat[, c('start','fu_days','event90','vlow', COV, AUTO), with=FALSE])
  meth <- rep('', ncol(X)); names(meth) <- names(X)
  for (v in c(COV, AUTO)) if (sum(is.na(X[[v]])) > 0) meth[v] <- 'pmm'
  sf <- survfit(Surv(start, fu_days, event90) ~ 1, data=X)
  X$nag <- approx(sf$time, sf$cumhaz, X$fu_days, rule=2, ties='ordered')$y
  X$nag[is.na(X$nag)] <- 0
  meth <- c(meth, nag = '')
  imp <- mice(X, m=20, maxit=10, method=meth, printFlag=FALSE, seed=20261105)
  fits <- lapply(seq_len(imp$m), function(k)
    coxph(Surv(start,fu_days,event90) ~ vlow + age + sex + sofa + charlson + lactate_filled2 + creat,
          data=complete(imp, k)))
  pl <- summary(pool(as.mira(fits)), conf.int=TRUE, exponentiate=TRUE)
  b <- pl[grep('^vlow', pl$term), ]
  cat(sprintf('%-46s n=%5d deaths=%4d | HR %.3f (%.3f-%.3f)\n', label, nrow(X), sum(X$event90),
              b$estimate, b[["2.5 %"]], b[["97.5 %"]]))
  data.frame(analysis=label, n=nrow(X), deaths=sum(X$event90), hr=b$estimate,
             lo=b[["2.5 %"]], hi=b[["97.5 %"]])
}

res <- list()
res[[1]] <- run(A, 'Primary (index TTE 0-72 h, first stay)')
res[[2]] <- run(d[alive], 'All in-window stays (first stay not required)')
res[[3]] <- run(A[echo_h <= 24], 'Index TTE within 24 h')
res[[4]] <- run(A[echo_h <= 48], 'Index TTE within 48 h')
q <- !(A$technical_quality %in% c('Poor','Suboptimal'))
cat(sprintf('  [quality] excluded %d of %d\n', sum(!q), nrow(A)))
res[[5]] <- run(A[q], 'Excluding poor or suboptimal image quality')
asx <- !grepl('^(Mod|Severe|Very severe|Low flow)', trimws(ifelse(is.na(A$aortic_stenosis), '', A$aortic_stenosis)))
cat(sprintf('  [aortic stenosis] excluded %d of %d\n', sum(!asx), nrow(A)))
res[[6]] <- run(A[asx], 'Excluding moderate-or-worse aortic stenosis')
pe <- !grepl('Moderate|Large', trimws(ifelse(is.na(A$pericardial_effusion), '', A$pericardial_effusion)))
cat(sprintf('  [effusion] excluded %d of %d\n', sum(!pe), nrow(A)))
res[[7]] <- run(A[pe], 'Excluding moderate or large pericardial effusion')

out <- do.call(rbind, res)
out$ci <- sprintf('%.2f (%.2f-%.2f)', out$hr, out$lo, out$hi)
write.csv(out, 'TableS4_mi.csv', row.names=FALSE)
cat('\nwritten TableS4_mi.csv\n'); print(out[, c('analysis','n','deaths','ci')], row.names=FALSE)
