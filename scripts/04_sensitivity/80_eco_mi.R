## 80_eco_mi.R — eCO (SV x heart rate, no BSA) instead of eCI, under multiple imputation
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(splines)); suppressMessages(library(mice))
suppressMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(x))
set.seed(20261105)
auc <- function(y, p) { r <- rank(p); n1 <- sum(y==1); n0 <- sum(y==0); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }

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
  !duplicated(stay_id), .(stay_id, measurement_id=as.character(measurement_id))],
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
A[, eco := ifelse(vti >= 5 & vti <= 45 & !is.na(vti_diam) & vti_diam >= 1.0 & vti_diam <= 3.2 &
                  sv >= 10 & sv <= 250 & !is.na(hrv), sv*hrv/1000, NA)]

AUX <- c('sapsii','mbp_min','hr_mean','temp_mean','spo2_min','bun','hgb','plt','wbc','bili',
         'inr','norepi','ventN','rrtN','vti_diam')
COV <- c('age','sex','sofa','charlson','lactate_filled2','creat')

H1 <- A[!is.na(lvef) & lvef > 0]
H2e <- H1[!is.na(eco)]
cat('H1', nrow(H1), 'deaths', sum(H1$event90), '| H2-eCO', nrow(H2e), 'deaths', sum(H2e$event90), '\n')

X <- as.data.frame(H2e[, c('start','fu_days','event90','vti','lvef','eco', COV, AUX), with=FALSE])
B <- ns(log(X$vti), df=3); colnames(B) <- paste0('v_s', 1:3); X <- cbind(X, B)
B <- ns(X$lvef, df=3); colnames(B) <- paste0('l_s', 1:3); X <- cbind(X, B)
B <- ns(log(X$eco), df=3); colnames(B) <- paste0('e_s', 1:3); X <- cbind(X, B)
meth <- rep('', ncol(X)); names(meth) <- names(X)
for (v in c(COV, AUX)) if (sum(is.na(X[[v]])) > 0) meth[v] <- 'pmm'
sf <- survfit(Surv(start, fu_days, event90) ~ 1, data=X)
X$nag <- approx(sf$time, sf$cumhaz, X$fu_days, rule=2, ties='ordered')$y
X$nag[is.na(X$nag)] <- 0
meth <- c(meth, nag='')
imp <- mice(X, m=20, maxit=10, method=meth, printFlag=FALSE, seed=20261105)
comp <- lapply(seq_len(imp$m), function(k) complete(imp, k))

forms <- list(clinical=paste('Surv(start,fu_days,event90) ~', paste(COV, collapse=' + ')),
              VTI_LVEF=paste('Surv(start,fu_days,event90) ~', paste(COV, collapse=' + '), '+ v_s1 + v_s2 + v_s3 + l_s1 + l_s2 + l_s3'),
              all=paste('Surv(start,fu_days,event90) ~', paste(COV, collapse=' + '),
                        '+ v_s1 + v_s2 + v_s3 + l_s1 + l_s2 + l_s3 + e_s1 + e_s2 + e_s3'))
P <- matrix(0, nrow(X), length(forms), dimnames=list(NULL, names(forms)))
for (nm in names(forms)) for (k in seq_along(comp)) {
  dk <- comp[[k]]
  f <- coxph(as.formula(forms[[nm]]), data=dk)
  bh <- basehaz(f, centered=FALSE)
  H90 <- approx(bh$time, bh$hazard, 90, rule=2)$y
  Hs <- approx(bh$time, bh$hazard, dk$start, rule=2)$y
  lp <- predict(f, newdata=dk, type='lp', reference='zero')
  P[, nm] <- P[, nm] + (1 - exp(-(H90 - Hs)*exp(lp)))/length(comp)
}
A0 <- sapply(colnames(P), function(v) auc(X$event90, P[, v]))
cat('\nAUC (eCO instead of eCI):\n'); print(round(A0, 4))
cat(sprintf('dAUC VTI+LVEF %+0.4f | dAUC all %+0.4f\n',
            A0['VTI_LVEF'] - A0['clinical'], A0['all'] - A0['clinical']))
write.csv(data.frame(model=names(A0), auc=as.numeric(A0),
                     dAUC=as.numeric(A0) - as.numeric(A0['clinical'])),
          'TableS5_eco_mi.csv', row.names=FALSE)
cat('written TableS5_eco_mi.csv\n')
