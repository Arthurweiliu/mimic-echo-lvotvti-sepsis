# 22_review4.R — 各切点 HR/CI + 验证期 HR
setwd('F:/openclaw_workspace/sepsis_echo_lit')
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA))
library(survival); cov<-'age+sex+sofa+charlson+lactate+creat'
dd<-d[complete.cases(d[,c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')])&d$fu_days>0&d$vti>0,]
cat('== 各切点 HR/CI (<cp vs >=cp), N=',nrow(dd),' ==\n')
for(cp in c(14,15,16,17,18,20)){
  dd$low<-as.numeric(dd$vti<cp)
  f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ low')),data=dd); b<-summary(f)$coefficients
  cat(sprintf('  <%d: HR=%.2f (%.2f-%.2f), n_low=%d\n', cp, exp(b['low',1]),exp(b['low',1]-1.96*b['low',3]),exp(b['low',1]+1.96*b['low',3]), sum(dd$low)))
}
cat('\n== 验证期(2014-2019)同模型 HR ==\n')
dd$vlow<-as.numeric(dd$vti<17)
va<-dd[dd$yg %in% c('2014 - 2016','2017 - 2019'),]
f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=va); b<-summary(f)$coefficients
cat(sprintf('  验证期 N=%d: VTI<17 HR=%.2f (%.2f-%.2f) p=%.3g\n', nrow(va), exp(b['vlow',1]),exp(b['vlow',1]-1.96*b['vlow',3]),exp(b['vlow',1]+1.96*b['vlow',3]),b['vlow',5]))
