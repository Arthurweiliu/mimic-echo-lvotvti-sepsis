# 19_selection.R — ③ TTE(±72h) vs 无 TTE 选择偏倚对照
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
cov<-read.csv('cohort_covariates.csv',stringsAsFactors=FALSE)
ak<-read.csv('aki_outcomes.csv',stringsAsFactors=FALSE)
lk<-read.csv('linked_cohort.csv',stringsAsFactors=FALSE)  # TTE ±72h stays
cov$tte<-as.numeric(cov$stay_id %in% lk$stay_id)
m<-merge(cov, ak[,c('stay_id','aki_any_7d','aki3_7d','crrt_7d')], by='stay_id', all.x=TRUE)
for(v in c('admission_age','sofa_fd','sapsii','charlson_comorbidity_index','lactate_max_fd','creatinine_max','norepi_eq_max_24h')) m[[v]]<-num(m[[v]])
T<-function(s){if(!s)return(NA);d<-as.Date(s);d}
m$int<-as.Date(m$icu_intime); m$dod<-as.Date(m$dod); m$disch<-as.Date(m$dischtime)
m$deathDate<-m$dod; m$deathDate[is.na(m$deathDate)&m$hospital_expire_flag==1]<-m$disch[is.na(m$deathDate)&m$hospital_expire_flag==1]
m$fu<-(m$deathDate-m$int)/1; m$ev90<-as.numeric(!is.na(m$fu)&m$fu<=90)
m$vent<-m$invasive_vent_24h %in% c('t','true',TRUE)
m$vaso<-as.numeric(!is.na(m$norepi_eq_max_24h) & m$norepi_eq_max_24h>0)
med<-function(x) sprintf('%.1f (%.1f-%.1f)', median(x,na.rm=TRUE), quantile(x,.25,na.rm=TRUE), quantile(x,.75,na.rm=TRUE))
pct<-function(x) sprintf('%d (%.0f%%)', sum(x,na.rm=TRUE), mean(x,na.rm=TRUE)*100)
pk<-function(v,g) tryCatch(format(kruskal.test(v~g)$p.value,digits=2,scientific=TRUE),error=function(e)'NA')
pc<-function(v,g) tryCatch(format(suppressWarnings(chisq.test(table(g,v)))$p.value,digits=2,scientific=TRUE),error=function(e)'NA')
grp<-m$tte
rows<-list()
r<-function(lbl,f,pf,fn_cont=TRUE) rows[[length(rows)+1]]<<-c(lbl, f(m[m$tte==1,]), f(m[m$tte==0,]), pf)
r('N', nrow, '', FALSE)
r('Age, years', function(x) med(x$admission_age), pk(m$admission_age,grp))
r('SOFA (first day)', function(x) med(x$sofa_fd), pk(m$sofa_fd,grp))
r('SAPS-II', function(x) med(x$sapsii), pk(m$sapsii,grp))
r('Charlson index', function(x) med(x$charlson_comorbidity_index), pk(m$charlson_comorbidity_index,grp))
r('Lactate, mmol/L', function(x) med(x$lactate_max_fd), pk(m$lactate_max_fd,grp))
r('Creatinine, mg/dL', function(x) med(x$creatinine_max), pk(m$creatinine_max,grp))
r('Invasive MV (24 h), n (%)', function(x) pct(x$vent), pc(m$vent,grp))
r('Vasopressor (24 h), n (%)', function(x) pct(x$vaso), pc(m$vaso,grp))
r('AKI (7 d), n (%)', function(x) pct(x$aki_any_7d=='t'), pc(m$aki_any_7d=='t',grp))
r('RRT (7 d), n (%)', function(x) pct(x$crrt_7d=='t'), pc(m$crrt_7d=='t',grp))
r('90-day mortality, n (%)', function(x) pct(x$ev90), pc(m$ev90,grp))
tab<-do.call(rbind,rows); colnames(tab)<-c('Variable','TTE \u00b172 h (n=14,136)','No TTE (n=18,835)','p')
write.csv(tab,'SELECTION_TTE_vs_noTTE.csv',row.names=FALSE)
cat('TTE:',sum(m$tte==1),' non-TTE:',sum(m$tte==0),'\n'); print(tab)
