# 37_audit_methods.R — methodological self-audit of the new analyses
{ .a <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x)); tr<-function(x) trimws(as.character(x))
d <- read.csv('h1_data.csv', stringsAsFactors=FALSE)
for (v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti','vti_diam','lvef','hr_mean')) d[[v]]<-num(d[[v]])
d$sex <- factor(ifelse(d$sex %in% c('M','F'), d$sex, NA)); d$stay_id<-as.character(d$stay_id)
suppressMessages(library(survival))
cov <- 'age+sex+sofa+charlson+lactate+creat'
ad <- read.csv('analysis_data.csv', stringsAsFactors=FALSE, colClasses=c(stay_id='character'))
ad <- ad[!duplicated(ad$stay_id), c('stay_id','measurement_id','body_surface_area','heart_rate_mean')]
d <- merge(d, ad, by='stay_id', all.x=TRUE); d$bsa<-num(d$body_surface_area); d$measurement_id<-as.character(d$measurement_id)
ex <- read.csv('echo_extra.csv', stringsAsFactors=FALSE); ex$measurement_id<-as.character(ex$measurement_id)
hrr <- ex[ex$measurement=='resting_hr', c('measurement_id','result')]; hrr<-hrr[!duplicated(hrr$measurement_id),]; names(hrr)[2]<-'rhr'
d <- merge(d, hrr, by='measurement_id', all.x=TRUE); d$rhr<-num(d$rhr)
af1<-read.csv('af_icu.csv',stringsAsFactors=FALSE); af2<-read.csv('af_icd.csv',stringsAsFactors=FALSE)
af1$stay_id<-as.character(af1$stay_id); af2$stay_id<-as.character(af2$stay_id)
d<-merge(d,af1[,c('stay_id','af_icu')],by='stay_id',all.x=TRUE); d<-merge(d,af2[,c('stay_id','af_icd')],by='stay_id',all.x=TRUE)
d$af<-as.numeric((!is.na(d$af_icu)&d$af_icu==1)|(!is.na(d$af_icd)&d$af_icd==1))

## ---- does the measurement table contain an independent CO/CI? ----
vd <- read.csv('variable_dictionary.csv', stringsAsFactors=FALSE)
hit <- grep('cardiac|output|index|\\bco\\b|stroke volume', vd$measurement_description, ignore.case=TRUE, value=TRUE)
cat('=== variable dictionary: any CO/CI/SV measurement? ===\n'); print(unique(hit)); cat('\n')

ccs <- c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')
A <- d[complete.cases(d[,ccs]) & num(d$fu_days)>0 & num(d$vti)>0, ]
A$sv <- A$vti*pi*(A$vti_diam/2)^2
A$hr <- ifelse(!is.na(A$rhr) & A$rhr>20 & A$rhr<220, A$rhr, NA)
A$co <- A$sv*A$hr/1000
A$ci <- ifelse(!is.na(A$bsa) & A$bsa>0.5, A$co/A$bsa, NA)
A$tachy <- as.numeric(A$hr>100)

## ---- (1) echo HR provenance ----
cat('=== (1) echo-time HR (resting_hr) sanity ===\n')
cat('available:', sum(!is.na(A$rhr)), sprintf('(%.1f%%)', 100*mean(!is.na(A$rhr))),
    '| median', round(median(A$rhr,na.rm=TRUE)), 'IQR', paste(round(quantile(A$rhr,c(.25,.75),na.rm=TRUE)),collapse='-'), '\n')
cat('Spearman(echo HR, ICU first-day mean HR):', round(cor(A$rhr, A$hr_mean, method='spearman', use='complete.obs'),2), '\n\n')

## ---- (2) LVOT diameter plausibility + sensitivity ----
cat('=== (2) LVOT diameter (cm) ===\n')
print(summary(A$vti_diam))
imp <- !is.na(A$vti_diam) & (A$vti_diam < 1.0 | A$vti_diam > 2.6)
cat('implausible (<1.0 or >2.6 cm):', sum(imp), sprintf('(%.1f%%)', 100*mean(imp)), '\n\n')

## ---- (3) analysis set & selection ----
Hok <- !is.na(A$sv)&A$sv>5&A$sv<250 & !is.na(A$co)&A$co>1&A$co<15 & !is.na(A$ci)&A$ci>0.5&A$ci<8
H <- A[Hok, ]; Ex <- A[!Hok, ]
cat('=== (3) haemodynamic subset vs excluded ===\n')
cat(sprintf('included n=%d (90d mort %.1f%%, SOFA %.1f, VTI %.1f) | excluded n=%d (mort %.1f%%, SOFA %.1f, VTI %.1f)\n',
  nrow(H), 100*mean(H$event90), median(H$sofa), median(H$vti), nrow(Ex), 100*mean(Ex$event90), median(Ex$sofa), median(Ex$vti)))
cat('diagnostic-model AUC for inclusion (age,sex,sofa,charlson,lactate,creat): ')
inc <- as.numeric(Hok); m <- glm(inc~age+sex+sofa+charlson+lactate+creat, data=A, family=binomial)
p<-predict(m,type='response'); r<-rank(p); n1<-sum(inc==1);n0<-sum(inc==0)
cat(round((sum(r[inc==1])-n1*(n1+1)/2)/(n1*n0),3), '\n\n')

## ---- (4) cancellation test: CO + HR joint ----
cat('=== (4) Is CO null because HR is inside it? ===\n')
g<-function(f,t,k=1){b<-summary(f)$coefficients; sprintf('%.3f (%.3f-%.3f) p=%.2g',exp(b[t,1]*k),exp((b[t,1]-1.96*b[t,3])*k),exp((b[t,1]+1.96*b[t,3])*k),b[t,5])}
cat('  CO alone        :', g(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ co')),data=H),'co'), '\n')
f2<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ co + hr')),data=H)
cat('  CO | HR adjusted: CO', g(f2,'co'), ' ; HR', g(f2,'hr',10), '\n')
cat('  CI alone        :', g(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ ci')),data=H),'ci'), '\n')
f3<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ ci + hr')),data=H)
cat('  CI | HR adjusted: CI', g(f3,'ci'), ' ; HR', g(f3,'hr',10), '\n\n')

## ---- (5) SV/CO associations under diameter restriction ----
cat('=== (5) sensitivity: restrict LVOT diameter to 1.0-2.6 cm ===\n')
Hr <- A[Hok & !imp, ]
for (v in c('vti','sv','co','ci')){ k <- c(vti=5,sv=10,co=1,ci=1)[v]
  cat(sprintf('  %-4s %s (n=%d)\n', v, g(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+',v)),data=Hr),v,k), nrow(Hr))) }
cat('\n')

## ---- (6) tachycardia interaction stability ----
cat('=== (6) VTI<17 x tachycardia across analysis sets ===\n')
for (nm in c('full primary (n=5410)','haemodynamic subset')){
  S <- if (grepl('full',nm)) A else H
  if (!'tachy' %in% names(S)) next
  S$vlow <- as.numeric(S$vti<17)
  f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow*tachy')),data=S); b<-summary(f)$coefficients
  it<-grep('^vlow:tachy$',rownames(b))
  cat(sprintf('  %-22s N=%5d  interaction %s\n', nm, nrow(S), g(f, rownames(b)[it])))
}
# AF scale-dependence
H$af <- A$af[match(H$stay_id, A$stay_id)]
cat('\n  AF x vlow(binary) :', g(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow*af')),data=transform(H,vlow=as.numeric(vti<17))),'vlow:af'), '\n')
H$x <- H$vti
cat('  AF x vti(continuous):', g(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ x*af')),data=H),'x:af'), '\n\n')

## ---- (7) C-index consistency ----
cat('=== (7) C-index (same subset, base vs +exposure) ===\n')
b0 <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov)), data=H)
cat('  base', round(summary(b0)$concordance[1],3), '| +VTI', round(summary(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti')),data=H))$concordance[1],3),
    '| +SV', round(summary(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ sv')),data=H))$concordance[1],3),
    '| +CO', round(summary(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ co')),data=H))$concordance[1],3), '\n')
cat('  (manuscript primary cohort C-index base 0.690 / +VTI 0.693)\n')
