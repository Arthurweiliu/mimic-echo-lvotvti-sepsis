## 42_baseline_smd.R — standardised mean differences for the two comparison tables
##   (S6) TTE ±72 h vs no TTE          (S8) LVOT-VTI available vs missing
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
smd_c <- function(v,g){ i<-!is.na(v)&!is.na(g); v<-v[i]; g<-g[i]
  m1<-mean(v[g==1]); m0<-mean(v[g==0]); s1<-sd(v[g==1]); s0<-sd(v[g==0])
  (m1-m0)/sqrt((s1^2+s0^2)/2) }
smd_b <- function(v,g){ i<-!is.na(v)&!is.na(g); v<-as.numeric(v[i]); g<-g[i]
  p1<-mean(v[g==1]); p0<-mean(v[g==0]); (p1-p0)/sqrt((p1*(1-p1)+p0*(1-p0))/2) }
med<-function(x) sprintf('%.1f (%.1f-%.1f)', median(x,na.rm=TRUE), quantile(x,.25,na.rm=TRUE), quantile(x,.75,na.rm=TRUE))
pct<-function(x) sprintf('%d (%.0f%%)', sum(x,na.rm=TRUE), mean(x,na.rm=TRUE)*100)

## ---------- S6: TTE vs no TTE ----------
cat('=== S6. TTE vs no TTE (SMD) ===\n')
cov<-read.csv('cohort_covariates.csv',stringsAsFactors=FALSE)
ak<-read.csv('aki_outcomes.csv',stringsAsFactors=FALSE)
lk<-read.csv('linked_cohort.csv',stringsAsFactors=FALSE)
cov$tte<-as.numeric(cov$stay_id %in% lk$stay_id)
m<-merge(cov, ak[,c('stay_id','aki_any_7d','crrt_7d')], by='stay_id', all.x=TRUE)
for(v in c('admission_age','sofa_fd','sapsii','charlson_comorbidity_index','lactate_max_fd','creatinine_max','norepi_eq_max_24h')) m[[v]]<-num(m[[v]])
m$int<-as.Date(m$icu_intime); m$dod<-as.Date(m$dod); m$disch<-as.Date(m$dischtime)
m$deathDate<-m$dod; m$deathDate[is.na(m$deathDate)&m$hospital_expire_flag==1]<-m$disch[is.na(m$deathDate)&m$hospital_expire_flag==1]
m$fu<-(m$deathDate-m$int); m$ev90<-as.numeric(!is.na(m$fu)&m$fu<=90)
m$vent<-m$invasive_vent_24h %in% c('t','true',TRUE); m$vaso<-as.numeric(!is.na(m$norepi_eq_max_24h)&m$norepi_eq_max_24h>0)
g<-m$tte
row6<-function(lbl,a,b,s) cat(sprintf('%-24s %-22s %-22s %s\n', lbl, a, b, s))
cat(sprintf('%-24s %-22s %-22s %s\n','Variable','TTE (n=14,136)','No TTE (n=18,835)','SMD'))
row6('N', as.character(sum(g==1)), as.character(sum(g==0)), '')
for(v in c('admission_age','sofa_fd','sapsii','charlson_comorbidity_index','lactate_max_fd','creatinine_max')){
  lab<-c(admission_age='Age, years',sofa_fd='SOFA (first day)',sapsii='SAPS-II',charlson_comorbidity_index='Charlson index',lactate_max_fd='Lactate, mmol/L',creatinine_max='Creatinine, mg/dL')[v]
  row6(lab, med(m[[v]][g==1]), med(m[[v]][g==0]), sprintf('%+.2f', smd_c(m[[v]],g))) }
row6('Invasive MV (24 h), n (%)', pct(m$vent[g==1]), pct(m$vent[g==0]), sprintf('%+.2f', smd_b(m$vent,g)))
row6('Vasopressor (24 h), n (%)', pct(m$vaso[g==1]), pct(m$vaso[g==0]), sprintf('%+.2f', smd_b(m$vaso,g)))
row6('AKI (7 d), n (%)', pct(m$aki_any_7d[g==1]=='t'), pct(m$aki_any_7d[g==0]=='t'), sprintf('%+.2f', smd_b(m$aki_any_7d=='t',g)))
row6('RRT (7 d), n (%)', pct(m$crrt_7d[g==1]=='t'), pct(m$crrt_7d[g==0]=='t'), sprintf('%+.2f', smd_b(m$crrt_7d=='t',g)))
row6('90-day mortality, n (%)', pct(m$ev90[g==1]), pct(m$ev90[g==0]), sprintf('%+.2f', smd_b(m$ev90,g)))

## ---------- S8: VTI available vs missing (first-ICU set, n = 12,120) ----------
cat('\n=== S8. LVOT-VTI available vs missing (SMD) ===\n')
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','sapsii','charlson','lactate','creat','norepi','event90','vti')) d[[v]]<-num(d[[v]])
d$ven<-d$vent %in% c('t','TRUE',TRUE); d$vas<-as.numeric(!is.na(d$norepi)&d$norepi>0)
av<-as.numeric(!is.na(d$vti)&d$vti>0)
cat(sprintf('%-24s %-22s %-22s %s\n','Variable','VTI available (n=7,253)','VTI missing (n=4,867)','SMD'))
cat(sprintf('%-24s %-22s %-22s %s\n','N', sum(av==1), sum(av==0), ''))
for(v in c('age','sofa','sapsii','charlson','lactate','creat')){
  lab<-c(age='Age, years',sofa='SOFA (first day)',sapsii='SAPS-II',charlson='Charlson index',lactate='Lactate, mmol/L',creat='Creatinine, mg/dL')[v]
  cat(sprintf('%-24s %-22s %-22s %s\n', lab, med(d[[v]][av==1]), med(d[[v]][av==0]), sprintf('%+.2f', smd_c(d[[v]],av)))) }
cat(sprintf('%-24s %-22s %-22s %s\n','Female, n (%)', sprintf('%d (%.0f%%)',sum(d$sex[av==1]=='F',na.rm=TRUE),100*mean(d$sex[av==1]=='F',na.rm=TRUE)),
  sprintf('%d (%.0f%%)',sum(d$sex[av==0]=='F',na.rm=TRUE),100*mean(d$sex[av==0]=='F',na.rm=TRUE)), sprintf('%+.2f', smd_b(d$sex=='F',av))))
cat(sprintf('%-24s %-22s %-22s %s\n','Invasive MV (24 h), n (%)', pct(d$ven[av==1]), pct(d$ven[av==0]), sprintf('%+.2f', smd_b(d$ven,av))))
cat(sprintf('%-24s %-22s %-22s %s\n','Vasopressor (24 h), n (%)', pct(d$vas[av==1]), pct(d$vas[av==0]), sprintf('%+.2f', smd_b(d$vas,av))))
cat(sprintf('%-24s %-22s %-22s %s\n','90-day mortality, n (%)', pct(d$event90[av==1]), pct(d$event90[av==0]), sprintf('%+.2f', smd_b(d$event90,av))))
