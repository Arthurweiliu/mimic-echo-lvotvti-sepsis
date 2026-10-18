## 41_review2_details.R — details requested by the second simulated review
##  (1) cut-point scan: n and HR for every scanned cut-point (fill the gaps in Table S1)
##  (2) IPTW: weight type, truncation, covariate balance (SMD before/after)
##  (3) the 27 deaths without a discharge/death timestamp
##  (4) direct LVOT-area sensitivity (is the Limitations claim supported?)
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
todt<-function(x){ x<-trimws(as.character(x)); x[x=='']<-NA; suppressWarnings(as.POSIXct(x, tz='UTC')) }
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','sapsii','charlson','lactate','creat','fu_days','event90','vti','vti_diam')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA)); d$stay_id<-as.character(d$stay_id)
suppressMessages(library(survival))
covs<-c('age','sex','sofa','charlson','lactate','creat'); cov<-paste(covs,collapse='+')
A<-d[complete.cases(d[,c('event90','fu_days',covs,'vti')]) & num(d$fu_days)>0 & num(d$vti)>0,]; A$vlow<-as.numeric(A$vti<17)

cat('=== (1) cut-point scan: n and HR for every cut-point ===\n')
for(cp in 13:26){
  x<-A; x$low<-as.numeric(x$vti<cp)
  m<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ low')),data=x); b<-summary(m)$coefficients
  cat(sprintf('  <%-2d cm  n_low=%4d  HR=%.2f (%.2f-%.2f) p=%.2g\n', cp, sum(x$low), exp(b['low',1]), exp(b['low',1]-1.96*b['low',3]), exp(b['low',1]+1.96*b['low',3]), b['low',5]))
}

cat('\n=== (2) IPTW weighting ===\n')
D<-d; D$inA<-as.numeric(D$stay_id %in% A$stay_id)
psc<-complete.cases(D[,c('age','sex','sofa','sapsii','charlson','lactate','creat')])
ps<-glm(inA ~ age+sex+sofa+sapsii+charlson+lactate+creat, data=D[psc,], family=binomial)
D$ps<-NA_real_; D$ps[psc]<-fitted(ps)
D$sw<-ifelse(D$inA==1, mean(D$inA[psc])/D$ps, (1-mean(D$inA[psc]))/(1-D$ps))   # stabilised
D$uw<-ifelse(D$inA==1, 1/D$ps, NA)                                            # unstabilised (analysis-set only)
cat('  unstabilised weights in Set A: range', paste(round(range(D$uw[D$inA==1],na.rm=TRUE),2),collapse=' - '),'\n')
cat('  stabilised weights in Set A  : range', paste(round(range(D$sw[D$inA==1],na.rm=TRUE),2),collapse=' - '),
    '| excluded:', paste(round(range(D$sw[D$inA==0],na.rm=TRUE),2),collapse=' - '),'\n')
smd<-function(v,grp,w=NULL){ i<-!is.na(v)&!is.na(grp); v<-v[i]; grp<-grp[i]
  if(is.null(w)) w<-rep(1,length(v)) else w<-w[i]
  m1<-weighted.mean(v[grp==1],w=w[grp==1]); m0<-weighted.mean(v[grp==0],w=w[grp==0])
  s1<-sqrt(sum(w[grp==1]*(v[grp==1]-m1)^2)/sum(w[grp==1])); s0<-sqrt(sum(w[grp==0]*(v[grp==0]-m0)^2)/sum(w[grp==0]))
  (m1-m0)/sqrt((s1^2+s0^2)/2) }
cat('  covariate balance (standardised mean differences, included vs excluded):\n')
base<-D[psc,]; cat(sprintf('    %-14s %8s %8s\n','covariate','unweighted','weighted'))
for(v in c('age','sofa','sapsii','charlson','lactate','creat')){
  cat(sprintf('    %-14s %8.3f %8.3f\n', v, smd(base[[v]], base$inA), smd(base[[v]], base$inA, w=base$sw))) }
cat(sprintf('    %-14s %8.3f %8.3f\n','female', smd(as.numeric(base$sex=='F'), base$inA), smd(as.numeric(base$sex=='F'), base$inA, w=base$sw)))
Aw<-A; Aw$sw<-D$sw[match(Aw$stay_id, D$stay_id)]
m<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=Aw,weights=sw,robust=TRUE); b<-summary(m)$coefficients
cat(sprintf('  stabilised IPTW HR=%.2f (%.2f-%.2f) p=%.2g\n', exp(b['vlow',1]), exp(b['vlow',1]-1.96*b['vlow',3]), exp(b['vlow',1]+1.96*b['vlow',3]), b['vlow',5]))

cat('\n=== (3) the 27 deaths with no discharge/death timestamp ===\n')
ad<-read.csv('analysis_data.csv',stringsAsFactors=FALSE,colClasses=c(stay_id='character'))
ad<-ad[!duplicated(ad$stay_id),c('stay_id','dod','dischtime','hospital_expire_flag')]
M<-merge(A[,c('stay_id','event90','hospital_death','fu_days')], ad, by='stay_id', all.x=TRUE)
dd<-M[M$event90==1,]
cat('  deaths in Set A:',nrow(dd),' | hospital_death=1:',sum(dd$hospital_death==1,na.rm=TRUE),
    ' | hospital_death=0:',sum(dd$hospital_death==0,na.rm=TRUE),' | hospital_death missing:',sum(is.na(dd$hospital_death)),'\n')
dd$dod2<-todt(dd$dod); dd$disch2<-todt(dd$dischtime)
cat('  with a usable dod AND discharge time:',sum(!is.na(dd$dod2)&!is.na(dd$disch2)),
    ' | hospital_death=1 among those:',sum(dd$hospital_death==1 & !is.na(dd$dod2)&!is.na(dd$disch2),na.rm=TRUE),'\n')
cat('  hospital_death flag distribution overall (Set A):\n'); print(table(A$hospital_death, useNA='ifany'))

cat('\n=== (4) direct LVOT-area sensitivity ===\n')
A2<-A; A2$diam_ok<-!is.na(A2$vti_diam) & A2$vti_diam>0
A2$sv<-A2$vti*pi*(A2$vti_diam/2)^2
g<-function(x,lab,frm,term='vlow'){ b<-summary(coxph(as.formula(frm),data=x))$coefficients
  cat(sprintf('  %-52s N=%5d  HR=%.2f (%.2f-%.2f) p=%.2g\n', lab, nrow(x), exp(b[term,1]), exp(b[term,1]-1.96*b[term,3]), exp(b[term,1]+1.96*b[term,3]), b[term,5])) }
g(A, 'VTI <17 vs >=17, all of Set A', paste('Surv(fu_days,event90)~',cov,'+ vlow'))
Ad<-A2[A2$diam_ok,]; Ad$vlow<-as.numeric(Ad$vti<17)
g(Ad, 'VTI <17 vs >=17, LVOT diameter recorded', paste('Surv(fu_days,event90)~',cov,'+ vlow'))
b<-summary(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ sv')),data=Ad))$coefficients
cat(sprintf('  %-52s N=%5d  HR=%.3f (%.3f-%.3f) p=%.2g  [per +10 mL]\n','derived stroke volume, LVOT diameter recorded',nrow(Ad),exp(b['sv',1]*10),exp((b['sv',1]-1.96*b['sv',3])*10),exp((b['sv',1]+1.96*b['sv',3])*10),b['sv',5]))
cat('\nDONE\n')
