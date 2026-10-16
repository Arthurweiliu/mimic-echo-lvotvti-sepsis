# 34_co_ci.R — does cardiac output (VTI x HR) beat VTI alone, and does the tachycardia
# interaction disappear once HR is incorporated?  Uses echo-time heart rate (resting_hr).
{ .a <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num <- function(x) suppressWarnings(as.numeric(x)); tr <- function(x) trimws(as.character(x))
d <- read.csv('h1_data.csv', stringsAsFactors = FALSE)
for (v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti','vti_diam','hr_mean')) d[[v]] <- num(d[[v]])
d$sex <- factor(ifelse(d$sex %in% c('M','F'), d$sex, NA))
d$stay_id <- as.character(d$stay_id)
suppressMessages(library(survival))
covs <- c('age','sex','sofa','charlson','lactate','creat'); cov <- paste(covs, collapse = '+')

## stay_id -> measurement_id  (h1_data has no measurement_id)
ad0 <- read.csv('analysis_data.csv', stringsAsFactors = FALSE, colClasses = c(stay_id = 'character'))
ad0 <- ad0[!duplicated(ad0$stay_id), c('stay_id','measurement_id','body_surface_area','heart_rate_mean')]
d <- merge(d, ad0, by = 'stay_id', all.x = TRUE); d$bsa <- num(d$body_surface_area)
d$measurement_id <- as.character(d$measurement_id)

## echo-time extras
ex <- read.csv('echo_extra.csv', stringsAsFactors = FALSE)
ex$measurement_id <- as.character(ex$measurement_id)
pick <- function(m, name) { s <- ex[ex$measurement == m, c('measurement_id','result')]
  s <- s[!duplicated(s$measurement_id), ]; names(s)[2] <- name; s }
d <- merge(d, pick('resting_hr','resting_hr'), by = 'measurement_id', all.x = TRUE)
d <- merge(d, pick('mitral_regurg','mr'), by = 'measurement_id', all.x = TRUE)
d <- merge(d, pick('tricuspid_regurg','tr_reg'), by = 'measurement_id', all.x = TRUE)
d$resting_hr <- num(d$resting_hr)

ccs <- c('event90','fu_days', covs, 'vti')
A <- d[complete.cases(d[, ccs]) & num(d$fu_days) > 0 & num(d$vti) > 0, ]
A$vlow <- as.numeric(A$vti < 17)
cat('Primary set N =', nrow(A), '| echo-time HR available:', sum(!is.na(A$resting_hr)),
    sprintf('(%.1f%%)', 100*mean(!is.na(A$resting_hr))), '\n')

## ---- derived haemodynamics ----
A$sv  <- A$vti * pi * (A$vti_diam/2)^2                       # mL
A$hr  <- ifelse(!is.na(A$resting_hr) & A$resting_hr > 20 & A$resting_hr < 220, A$resting_hr, NA)
A$co  <- A$sv * A$hr / 1000                                  # L/min
A$ci  <- ifelse(!is.na(A$bsa) & A$bsa > 0.5, A$co / A$bsa, NA)  # L/min/m2
ok <- !is.na(A$sv) & A$sv>5 & A$sv<250 & !is.na(A$co) & A$co>1 & A$co<15 & !is.na(A$ci) & A$ci>0.5 & A$ci<8
H <- A[ok, ]
cat('Haemodynamic-complete subset:', nrow(H), sprintf('(%.1f%%)', 100*nrow(H)/nrow(A)), '\n')
cat('  echo HR med', round(median(H$hr)), '| SV med', round(median(H$sv)), 'mL | CO med',
    round(median(H$co),2), 'L/min | CI med', round(median(H$ci),2), 'L/min/m2\n')
cat('  Spearman: VTI-CO', round(cor(H$vti,H$co,method='spearman'),2),
    '| VTI-HR', round(cor(H$vti,H$hr,method='spearman'),2),
    '| SV-HR', round(cor(H$sv,H$hr,method='spearman'),2), '\n\n')

bb <- function(f, t, k) { b <- summary(f)$coefficients
  sprintf('%.3f (%.3f-%.3f) p=%.2g', exp(b[t,1]*k), exp((b[t,1]-1.96*b[t,3])*k), exp((b[t,1]+1.96*b[t,3])*k), b[t,5]) }
cat('=== Adjusted association with 90-day mortality (haemodynamic subset) ===\n')
fits <- list(
  'VTI  per +5 cm'      = coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti')), data=H),
  'SV   per +10 mL'     = coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ sv')),  data=H),
  'CO   per +1 L/min'   = coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ co')),  data=H),
  'CI   per +1 L/min/m2'= coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ ci')),  data=H))
for (nm in names(fits)) cat(sprintf('  %-20s %s\n', nm, bb(fits[[nm]], c('vti','sv','co','ci')[match(nm,names(fits))], 1)))
cat('  C-index:'); for (nm in names(fits)) cat(sprintf('  %s=%.3f', sub(' .*','',nm), summary(fits[[nm]])$concordance[1])); cat('\n')
H$ci_low <- as.numeric(H$ci < 2.2)
f <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ ci_low')), data=H); b<-summary(f)$coefficients
cat(sprintf('  Low CI (<2.2 L/min/m2), n=%d: HR=%.2f (%.2f-%.2f) p=%.2g\n', sum(H$ci_low),
  exp(b['ci_low',1]), exp(b['ci_low',1]-1.96*b['ci_low',3]), exp(b['ci_low',1]+1.96*b['ci_low',3]), b['ci_low',5]))

cat('\n=== Does the tachycardia interaction survive when HR is built INTO the exposure? ===\n')
H$tachy <- as.numeric(H$hr > 100)
g <- function(f, t) { b <- summary(f)$coefficients
  sprintf('%.3f (%.3f-%.3f) p=%.3g', exp(b[t,1]), exp(b[t,1]-1.96*b[t,3]), exp(b[t,1]+1.96*b[t,3]), b[t,5]) }
f1 <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow*tachy')), data=H)
cat('  (i) VTI<17  x tachycardia :', g(f1,'vlow:tachy'), '\n')
f2 <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ co*tachy')), data=H)
cat('  (ii) CO x tachycardia     :', g(f2,'co:tachy'), '\n')
f3 <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ ci*tachy')), data=H)
cat('  (iii) CI x tachycardia    :', g(f3,'ci:tachy'), '\n')
cat('  CO by stratum: HR<=100 ', bb(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ co')), data=H[H$tachy==0,]),'co',1),
    ' | HR>100 ', bb(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ co')), data=H[H$tachy==1,]),'co',1), '\n')

cat('\n=== MR / TR added to the modifier panel (90-day) ===\n')
A2 <- A
A2$mr_mod <- grepl('mod', tr(A2$mr), ignore.case = TRUE) & !grepl('mild', tr(A2$mr), ignore.case = TRUE)
A2$tr_mod <- grepl('mod|severe', tr(A2$tr_reg), ignore.case = TRUE) & !grepl('mild|trace', tr(A2$tr_reg), ignore.case = TRUE)
A2$mr_mod[is.na(A2$mr_mod)] <- FALSE; A2$tr_mod[is.na(A2$tr_mod)] <- FALSE
for (m in c('mr_mod','tr_mod')){
  s <- A2[complete.cases(A2[, ccs]), ]
  f <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow*',m)), data=s); b<-summary(f)$coefficients
  it <- grep(paste0('^vlow:',m,'TRUE$'), rownames(b))
  cat(sprintf('  %-7s preval %.1f%% | interaction ', m, 100*mean(s[[m]])), g(f, rownames(b)[it]), '\n')
}
cat('\nDONE\n')
