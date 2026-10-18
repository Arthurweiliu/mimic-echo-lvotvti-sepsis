## 40_reviewer_sensitivity.R — analyses requested by the simulated review
##  (a) selection into the analysis set: baseline + inverse-probability-of-selection weighting
##  (b) measurement provenance: LVOT diameter availability; echo-time HR comes from the same study
##  (c) outcome ascertainment: how much of the 90-day mortality is post-discharge?
##  (d) time-segmented hazard ratios (0-7 / 7-28 / 28-90 days)
##  (e) derived-CO error simulation
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
todt<-function(x){ x<-trimws(as.character(x)); x[x=='']<-NA; suppressWarnings(as.POSIXct(x, tz='UTC')) }
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','sapsii','charlson','lactate','creat','fu_days','event90','vti','vti_diam','lvef')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA)); d$stay_id<-as.character(d$stay_id)
suppressMessages(library(survival))
covs<-c('age','sex','sofa','charlson','lactate','creat'); cov<-paste(covs,collapse='+')
ccs<-c('event90','fu_days',covs,'vti')
A<-d[complete.cases(d[,ccs]) & num(d$fu_days)>0 & num(d$vti)>0,]; A$vlow<-as.numeric(A$vti<17)
cat('Set A N =',nrow(A),' of',nrow(d),'patients with VTI available\n\n')

hr <- function(f,lab,term='vlow'){ b<-summary(f)$coefficients
  cat(sprintf('%-46s HR=%.2f (%.2f-%.2f) p=%.2g\n', lab, exp(b[term,1]), exp(b[term,1]-1.96*b[term,3]), exp(b[term,1]+1.96*b[term,3]), b[term,5])) }

## (a) selection: Set A vs excluded (among the 12,120 first-ICU patients)
D<-d; D$inA<-as.numeric(D$stay_id %in% A$stay_id)
psc<-complete.cases(D[,c('age','sex','sofa','sapsii','charlson','lactate','creat')])
ps<-glm(inA ~ age+sex+sofa+sapsii+charlson+lactate+creat, data=D[psc,], family=binomial)
D$ps<-NA_real_; D$ps[psc]<-fitted(ps)
D$w<-ifelse(D$inA==1, 1/D$ps, NA)
cat('=== (a) selection into Set A ===\n')
cat('  mean propensity score in Set A =',round(mean(D$ps[D$inA==1]),3),
    '| in excluded =',round(mean(D$ps[D$inA==0],na.rm=TRUE),3),'\n')
cat('  weight range (Set A):',paste(round(range(D$w[D$inA==1]),2),collapse=' - '),'\n')
Aw<-A; Aw$w<-D$w[match(Aw$stay_id, D$stay_id)]
hr(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=Aw,weights=w,robust=TRUE), 'IPTW (weighted to the 12,120 cohort)')
hr(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=Aw[Aw$w<=10,],weights=w,robust=TRUE), 'IPTW, weights truncated at 10')
hr(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=A), 'unweighted (reference)')

## (b) measurement provenance
cat('\n=== (b) measurement provenance ===\n')
cat('  LVOT diameter recorded in Set A:',sum(!is.na(A$vti_diam)),sprintf('(%.1f%%)',100*mean(!is.na(A$vti_diam))),'\n')
ex<-read.csv('echo_extra.csv',stringsAsFactors=FALSE)
cat('  echo_extra measurements available:',paste(setdiff(sort(unique(ex$measurement)),''),collapse=', '),'\n')
cat('  resting_hr is keyed by measurement_id (same echo study as lvot_vti):', 'resting_hr' %in% ex$measurement, '\n')

## (c) outcome ascertainment: in-hospital vs post-discharge deaths
cat('\n=== (c) where do the 90-day deaths occur? ===\n')
ad<-read.csv('analysis_data.csv',stringsAsFactors=FALSE,colClasses=c(stay_id='character'))
ad<-ad[!duplicated(ad$stay_id),c('stay_id','dod','dischtime','icu_intime','hospital_expire_flag')]
ad$dod<-todt(ad$dod); ad$disch<-todt(ad$dischtime); ad$int<-todt(ad$icu_intime)
ad$post<-as.numeric(!is.na(ad$dod) & !is.na(ad$disch) & ad$dod>ad$disch)
M<-merge(A[,c('stay_id','event90','hospital_death')], ad, by='stay_id', all.x=TRUE)
dd<-M[M$event90==1,]
cat(sprintf('  deaths within 90 days: %d | in-hospital: %d (%.1f%%) | after discharge: %d (%.1f%%)\n',
  nrow(dd), sum(dd$hospital_death==1,na.rm=TRUE), 100*mean(dd$hospital_death==1,na.rm=TRUE),
  sum(dd$post,na.rm=TRUE), 100*mean(dd$post,na.rm=TRUE)))
hr(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=A), '  90-day mortality (primary)')
Aih<-A; Aih$ev_ih<-as.numeric(Aih$hospital_death==1)
hr(coxph(as.formula(paste('Surv(pmin(fu_days,90),ev_ih)~',cov,'+ vlow')),data=Aih), '  in-hospital mortality as the outcome')

## (d) time-segmented hazard ratios
cat('\n=== (d) time-segmented HR (period-specific) ===\n')
sp<-survSplit(Surv(fu_days,event90)~., data=A, cut=c(7,28), episode='period', end='tstop', start='tstart')
sp$period<-factor(sp$period, labels=c('0-7 d','7-28 d','28-90 d'))
f2<-coxph(as.formula(paste('Surv(tstart,tstop,event90)~',cov,'+ vlow*period')),data=sp)
b<-summary(f2)$coefficients
for(i in rownames(b)) if(grepl('vlow',i)) cat(sprintf('  %-28s HR=%.2f (%.2f-%.2f) p=%.2g\n', i, exp(b[i,1]), exp(b[i,1]-1.96*b[i,3]), exp(b[i,1]+1.96*b[i,3]), b[i,5]))

## (e) derived-CO error simulation
cat('\n=== (e) CO derivation error simulation (Set C) ===\n')
ad0<-read.csv('analysis_data.csv',stringsAsFactors=FALSE,colClasses=c(stay_id='character'))
ad0<-ad0[!duplicated(ad0$stay_id),c('stay_id','measurement_id','body_surface_area')]
A2<-merge(A, ad0, by='stay_id', all.x=TRUE)
A2$bsa<-num(A2$body_surface_area); A2$measurement_id<-as.character(A2$measurement_id)
s<-ex[ex$measurement=='resting_hr',c('measurement_id','result')]; s<-s[!duplicated(s$measurement_id),]; names(s)[2]<-'resting_hr'
A2<-merge(A2,s,by='measurement_id',all.x=TRUE); A2$resting_hr<-num(A2$resting_hr)
A2$sv<-A2$vti*pi*(A2$vti_diam/2)^2; A2$hr<-ifelse(A2$resting_hr>20 & A2$resting_hr<220, A2$resting_hr, NA)
A2$co<-A2$sv*A2$hr/1000; A2$ci<-ifelse(!is.na(A2$bsa)&A2$bsa>0.5, A2$co/A2$bsa, NA)
H<-A2[!is.na(A2$sv)&A2$sv>5&A2$sv<250&!is.na(A2$co)&A2$co>1&A2$co<15&!is.na(A2$ci)&A2$ci>0.5&A2$ci<8,]
cat('  Set C N =',nrow(H),'\n')
set.seed(42)
res<-t(replicate(200, {
  vti_e <- H$vti * rnorm(nrow(H), 1, 0.05)
  diam_e<- H$vti_diam * rnorm(nrow(H), 1, 0.05)
  H$co_e<- vti_e*pi*(diam_e/2)^2 * H$hr/1000
  b<-summary(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ co_e')),data=H))$coefficients
  c(est=exp(b['co_e',1]), lo=exp(b['co_e',1]-1.96*b['co_e',3]), hi=exp(b['co_e',1]+1.96*b['co_e',3]))
}))
cat(sprintf('  CO per +1 L/min with 5%% error in VTI and LVOT diameter: mean HR=%.3f (mean CI %.3f-%.3f)\n',
  mean(res[,'est']), mean(res[,'lo']), mean(res[,'hi'])))
b0<-summary(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ co')),data=H))$coefficients
cat(sprintf('  observed CO per +1 L/min: HR=%.3f (%.3f-%.3f) p=%.2g\n', exp(b0['co',1]), exp(b0['co',1]-1.96*b0['co',3]), exp(b0['co',1]+1.96*b0['co',3]), b0['co',5]))
cat('\nDONE\n')
