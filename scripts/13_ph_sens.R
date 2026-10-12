# 13_ph_sens.R — PH 假设违背的敏感性（时间变换 tt + 限制随访）
setwd('F:/openclaw_workspace/sepsis_echo_lit')
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA))
suppressMessages(library(survival))
dd<-d[complete.cases(d[,c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')])&d$fu_days>0&d$vti>0,]
dd$vlow<-as.numeric(dd$vti<17)
cat('N=',nrow(dd),'\n')
# 1) 限制随访 28 天
dd$e28<-as.numeric(dd$fu_days<=28 & dd$event90==1); dd$t28<-pmin(dd$fu_days,28)
f2<-coxph(Surv(t28,e28)~age+sex+sofa+charlson+lactate+creat+vlow, data=dd); b2<-summary(f2)$coefficients
cat(sprintf('限制随访28天: VTI<17 HR=%.2f (%.2f-%.2f) p=%.3g\n', exp(b2['vlow',1]), exp(b2['vlow',1]-1.96*b2['vlow',3]), exp(b2['vlow',1]+1.96*b2['vlow',3]), b2['vlow',5]))
# 2) 时间变换（VTI<17 效应随时间变化）
f3<-coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+vlow+tt(vlow),
          tt=function(x,t,...) x*log(t+1), data=dd)
b3<-summary(f3)$coefficients
cat(sprintf('tt 模型: vlow 基线项 HR=%.2f (%.2f-%.2f)\n', exp(b3['vlow',1]), exp(b3['vlow',1]-1.96*b3['vlow',3]), exp(b3['vlow',1]+1.96*b3['vlow',3])))
cat(sprintf('        vlow:log(t+1) 系数=%.4f (p=%.3g) → 效应随时间%s\n', b3['tt(vlow)',1], b3['tt(vlow)',5], ifelse(b3['tt(vlow)',5]<0.05,'显著变化','基本稳定')))
