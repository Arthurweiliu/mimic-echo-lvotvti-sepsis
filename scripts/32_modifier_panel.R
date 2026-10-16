# 32_modifier_panel.R — systematic "validity / reliability modifier" panel for LVOT-VTI
# Tests whether the prognostic association (vlow <17 vs >=17) is modified by:
#   rhythm (AF), valve disease (AS>=mod, AR>=mod), tachycardia, image quality
# and re-tests the exposure as AREA-CORRECTED stroke volume (condition B).
{ .a <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num <- function(x) suppressWarnings(as.numeric(x)); tr <- function(x) trimws(as.character(x))
d <- read.csv('h1_data.csv', stringsAsFactors = FALSE)
for (v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti','vti_diam','hr_mean'))
  d[[v]] <- num(d[[v]])
d$sex <- factor(ifelse(d$sex %in% c('M','F'), d$sex, NA))
suppressMessages(library(survival))
covs <- c('age','sex','sofa','charlson','lactate','creat'); cov <- paste(covs, collapse = '+')

## AF flag
af1 <- read.csv('af_icu.csv', stringsAsFactors = FALSE); af2 <- read.csv('af_icd.csv', stringsAsFactors = FALSE)
d$stay_id <- as.character(d$stay_id); af1$stay_id <- as.character(af1$stay_id); af2$stay_id <- as.character(af2$stay_id)
d <- merge(d, af1[, c('stay_id','af_icu')], by = 'stay_id', all.x = TRUE)
d <- merge(d, af2[, c('stay_id','af_icd')], by = 'stay_id', all.x = TRUE)
d$af <- as.numeric((!is.na(d$af_icu) & d$af_icu == 1) | (!is.na(d$af_icd) & d$af_icd == 1))

## BSA (from analysis_data) for indexing SV
ad <- read.csv('analysis_data.csv', stringsAsFactors = FALSE, colClasses = c(stay_id = 'character'))
ad <- ad[!duplicated(ad$stay_id), c('stay_id','body_surface_area')]
d <- merge(d, ad, by = 'stay_id', all.x = TRUE); d$bsa <- num(d$body_surface_area)

## modifiers
d$as_mod  <- grepl('^mod|mod to severe|^severe|very severe|low flow', tr(d$aortic_stenosis), ignore.case = TRUE)
d$ar_mod  <- grepl('mod \\[2|mod-severe|severe \\[4', tr(d$aortic_regurg), ignore.case = TRUE)
d$tech_sub<- tr(d$technical_quality) %in% c('Suboptimal','Poor')
d$tachy   <- !is.na(d$hr_mean) & d$hr_mean > 100
for (v in c('as_mod','ar_mod','tech_sub','tachy')) d[[v]][is.na(d[[v]])] <- FALSE

ccs <- c('event90','fu_days', covs, 'vti')
A <- d[complete.cases(d[, ccs]) & num(d$fu_days) > 0 & num(d$vti) > 0, ]
A$vlow <- as.numeric(A$vti < 17)
cat('Primary set N =', nrow(A), ' events =', sum(A$event90 == 1), '\n')
cat(sprintf('Modifier prevalence: AF %.1f%% | AS>=mod %.1f%% | AR>=mod %.1f%% | tachy(>100) %.1f%% | tech sub/poor %.1f%%\n\n',
  100*mean(A$af), 100*mean(A$as_mod), 100*mean(A$ar_mod), 100*mean(A$tachy), 100*mean(A$tech_sub)))

## interaction panel (binary vlow)
panel <- function(mod, lab) {
  f <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ vlow*', mod)), data = A)
  b <- summary(f)$coefficients; it <- grep(paste0('^vlow:', mod, 'TRUE$|^vlow:', mod, '$'), rownames(b))
  it <- it[1]
  hr0 <- summary(coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ vlow')), data = A[A[[mod]] == 0, ]))$coefficients
  hr1 <- summary(coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ vlow')), data = A[A[[mod]] == 1, ]))$coefficients
  data.frame(Modifier = lab, n_absent = sum(A[[mod]] == 0), n_present = sum(A[[mod]] == 1),
    HR_absent = round(exp(hr0['vlow',1]),2), abs_lo = round(exp(hr0['vlow',1]-1.96*hr0['vlow',3]),2), abs_hi = round(exp(hr0['vlow',1]+1.96*hr0['vlow',3]),2),
    HR_present = round(exp(hr1['vlow',1]),2), pre_lo = round(exp(hr1['vlow',1]-1.96*hr1['vlow',3]),2), pre_hi = round(exp(hr1['vlow',1]+1.96*hr1['vlow',3]),2),
    int_HR = round(exp(b[it,1]),2), int_lo = round(exp(b[it,1]-1.96*b[it,3]),2), int_hi = round(exp(b[it,1]+1.96*b[it,3]),2), int_p = b[it,5])
}
res <- do.call(rbind, list(
  panel('af','Atrial fibrillation'), panel('as_mod','Aortic stenosis >= mod'),
  panel('ar_mod','Aortic regurgitation >= mod'), panel('tachy','Tachycardia HR>100'),
  panel('tech_sub','Image quality suboptimal/poor')))
res$int_p_FDR <- round(p.adjust(res$int_p, method = 'BH'), 3); res$int_p <- signif(res$int_p, 2)
cat('=== Effect-modification panel (vlow <17 vs >=17, 90-day mortality) ===\n'); print(res, row.names = FALSE)
cat('\n(prespecified in manuscript: AF, AS>=mod, image quality; AR and tachycardia are new)\n')
write.csv(res, 'Item3_modifier_panel.csv', row.names = FALSE)

## ---- condition (B): area-corrected stroke volume ----
cat('\n=== Condition B: LVOT area-corrected stroke volume ===\n')
cat('vti_diam range:', paste(range(A$vti_diam, na.rm = TRUE), collapse = '-'), 'cm ; non-missing', sum(!is.na(A$vti_diam)), '\n')
A$sv <- A$vti * pi * (A$vti_diam/2)^2       # mL
A$svi <- A$sv / A$bsa                       # mL/m2
ok <- !is.na(A$sv) & A$sv > 5 & A$sv < 250
cat('SV computable & plausible:', sum(ok), sprintf('(%.1f%%)', 100*mean(ok)), '; median SV', round(median(A$sv[ok])), 'mL [IQR',
    paste(round(quantile(A$sv[ok], c(.25,.75))), collapse='-'), '] ; median SVI', round(median(A$svi[ok], na.rm=TRUE)), 'mL/m2\n')
Ak <- A[ok, ]
f_vti <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ vti')), data = Ak)
f_sv  <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ sv')),  data = Ak)
f_svi <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ svi')), data = Ak)
bb <- function(f, t, k) { b <- summary(f)$coefficients
  sprintf('%.3f (%.3f-%.3f) p=%.2g', exp(b[t,1]*k), exp((b[t,1]-1.96*b[t,3])*k), exp((b[t,1]+1.96*b[t,3])*k), b[t,5]) }
cat('  per +5 cm VTI      :', bb(f_vti,'vti',5), '\n')
cat('  per +10 mL SV      :', bb(f_sv,'sv',10), '\n')
cat('  per +5 mL/m2 SVI   :', bb(f_svi,'svi',5), '\n')
cat(sprintf('  C-index: VTI %.3f | SV %.3f | SVI %.3f\n', summary(f_vti)$concordance[1], summary(f_sv)$concordance[1], summary(f_svi)$concordance[1]))
cat('  LVOT diam vs VTI correlation (Spearman):', round(cor(Ak$vti_diam, Ak$vti, method='spearman', use='complete.obs'), 3), '\n')
## AF x (vti vs sv) interaction
for (ex in c('vti','sv')){
  Ak$ex <- Ak[[ex]]
  f <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ ex*af')), data = Ak); b <- summary(f)$coefficients
  it <- grep('^ex:af$', rownames(b))
  cat(sprintf('  Interaction AF x %-4s : HR=%.3f (%.3f-%.3f) p=%.3f\n', ex, exp(b[it,1]), exp(b[it,1]-1.96*b[it,3]), exp(b[it,1]+1.96*b[it,3]), b[it,5]))
}
