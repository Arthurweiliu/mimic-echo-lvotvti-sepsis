# 33_modifier_joint.R — corrected condition-B + joint multi-modifier analysis
{ .a <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num <- function(x) suppressWarnings(as.numeric(x)); tr <- function(x) trimws(as.character(x))
d <- read.csv('h1_data.csv', stringsAsFactors = FALSE)
for (v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti','vti_diam','hr_mean')) d[[v]] <- num(d[[v]])
d$sex <- factor(ifelse(d$sex %in% c('M','F'), d$sex, NA))
suppressMessages(library(survival))
covs <- c('age','sex','sofa','charlson','lactate','creat'); cov <- paste(covs, collapse = '+')
af1 <- read.csv('af_icu.csv', stringsAsFactors = FALSE); af2 <- read.csv('af_icd.csv', stringsAsFactors = FALSE)
d$stay_id <- as.character(d$stay_id); af1$stay_id <- as.character(af1$stay_id); af2$stay_id <- as.character(af2$stay_id)
d <- merge(d, af1[, c('stay_id','af_icu')], by = 'stay_id', all.x = TRUE)
d <- merge(d, af2[, c('stay_id','af_icd')], by = 'stay_id', all.x = TRUE)
d$af <- as.numeric((!is.na(d$af_icu) & d$af_icu == 1) | (!is.na(d$af_icd) & d$af_icd == 1))
ad <- read.csv('analysis_data.csv', stringsAsFactors = FALSE, colClasses = c(stay_id = 'character'))
ad <- ad[!duplicated(ad$stay_id), c('stay_id','body_surface_area')]
d <- merge(d, ad, by = 'stay_id', all.x = TRUE); d$bsa <- num(d$body_surface_area)
d$tachy <- !is.na(d$hr_mean) & d$hr_mean > 100; d$tachy[is.na(d$tachy)] <- FALSE
ccs <- c('event90','fu_days', covs, 'vti')
A <- d[complete.cases(d[, ccs]) & num(d$fu_days) > 0 & num(d$vti) > 0, ]
A$vlow <- as.numeric(A$vti < 17)

hr_get <- function(f, t) { b <- summary(f)$coefficients
  sprintf('%.2f (%.2f-%.2f) p=%.3g', exp(b[t,1]), exp(b[t,1]-1.96*b[t,3]), exp(b[t,1]+1.96*b[t,3]), b[t,5]) }

cat('=== (A) Joint AF + tachycardia interaction (mutually adjusted) ===\n')
f <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ vlow*af + vlow*tachy')), data = A)
cat('  vlow x AF    :', hr_get(f, 'vlow:af'), '\n')
cat('  vlow x tachy :', hr_get(f, 'vlow:tachyTRUE'), '\n')

cat('\n=== (B) AF effect within HR strata ===\n')
for (g in c('HR<=100','HR>100')) {
  s <- if (g == 'HR<=100') A[!A$tachy, ] else A[A$tachy, ]
  if (nrow(s) > 200) { f2 <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ vlow*af')), data = s)
    cat(sprintf('  %-8s N=%5d : vlow x AF %s\n', g, nrow(s), hr_get(f2, 'vlow:af'))) }
}

cat('\n=== (C) Continuous HR as modifier ===\n')
f3 <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ vlow*hr_mean')), data = A)
cat('  vlow x hr_mean(per 10 bpm):', hr_get(f3, 'vlow:hr_mean'), '\n')

cat('\n=== (D) Condition B: area-corrected SV / SVI ===\n')
A$sv  <- A$vti * pi * (A$vti_diam/2)^2
A$svi <- ifelse(!is.na(A$bsa) & A$bsa > 0.5, A$sv / A$bsa, NA)
ok <- !is.na(A$sv) & A$sv > 5 & A$sv < 250 & !is.na(A$svi) & A$svi > 5 & A$svi < 120
Ak <- A[ok, ]
cat('  SV & SVI computable:', nrow(Ak), sprintf('(%.1f%% of primary set)', 100*nrow(Ak)/nrow(A)),
    '| median SVI', round(median(Ak$svi)), 'mL/m2\n')
bb <- function(f, t, k) { b <- summary(f)$coefficients
  sprintf('%.3f (%.3f-%.3f) p=%.2g', exp(b[t,1]*k), exp((b[t,1]-1.96*b[t,3])*k), exp((b[t,1]+1.96*b[t,3])*k), b[t,5]) }
f_vti <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ vti')), data = Ak)
f_sv  <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ sv')),  data = Ak)
f_svi <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ svi')), data = Ak)
cat('  per +5 cm VTI    :', bb(f_vti,'vti',5), '\n')
cat('  per +10 mL SV    :', bb(f_sv,'sv',10), '\n')
cat('  per +5 mL/m2 SVI :', bb(f_svi,'svi',5), '\n')
cat(sprintf('  C-index: VTI %.3f | SV %.3f | SVI %.3f\n',
    summary(f_vti)$concordance[1], summary(f_sv)$concordance[1], summary(f_svi)$concordance[1]))
cat('  Spearman(VTI, LVOT diam) =', round(cor(Ak$vti, Ak$vti_diam, method='spearman', use='complete.obs'),3),
    '| Spearman(VTI, SV) =', round(cor(Ak$vti, Ak$sv, method='spearman', use='complete.obs'),3), '\n')
for (ex in c('vti','sv')){
  Ak$ex <- Ak[[ex]]
  f4 <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ ex*af')), data = Ak)
  cat(sprintf('  AF x %-4s interaction:', ex), hr_get(f4, 'ex:af'), '\n')
}
cat('\nDONE\n')
