# 17_reviewer_checks.R — 针对审稿意见的补充计算
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)  # first-ICU
for(v in c('age','sofa','sapsii','charlson','lactate','creat','fu_days','event90','death28','vti')) d[[v]]<-num(d[[v]])
library(survival)
cat('== 缺失比例 (first-ICU 队列 N=',nrow(d),') ==\n')
for(v in c('vti','lactate','creat','sofa','charlson','age','sex')) cat(sprintf('  %-8s missing %d (%.1f%%)\n', v, sum(is.na(d[[v]])|(v=='vti'&d[[v]]==0)), mean(is.na(d[[v]])|(v=='vti'&d[[v]]==0))*100))
cat('  VTI 缺失 in TTE-eligible, first-ICU:', sum(is.na(d$vti)|d$vti==0), 'of', nrow(d), '\n')
cat('\n== 28 天死亡 (first-ICU) ==\n')
d$vlow<-as.numeric(d$vti<17)
dd<-d[complete.cases(d[,c('death28','fu_days','age','sex','sofa','charlson','lactate','creat','vti')])&d$fu_days>0&d$vti>0,]
dd$t28<-pmin(dd$fu_days,28); dd$e28<-as.numeric(dd$fu_days<=28 & dd$event90==1)
f<-coxph(Surv(t28,e28)~age+sex+sofa+charlson+lactate+creat+vlow,data=dd); b<-summary(f)$coefficients
cat(sprintf('  VTI<17 vs >=17: HR=%.2f (%.2f-%.2f) p=%.3g | 28d 死亡率 <17=%.1f%% vs >=17=%.1f%%\n',
  exp(b['vlow',1]),exp(b['vlow',1]-1.96*b['vlow',3]),exp(b['vlow',1]+1.96*b['vlow',3]),b['vlow',5],
  mean(dd$e28[dd$vlow==1])*100, mean(dd$e28[dd$vlow==0])*100))
cat('\n== Landmark 分析 (自 72h 起算) ==\n')
lm<-dd[dd$fu_days>3,]; lm$t2<-lm$fu_days-3; lm$e2<-lm$event90
f2<-coxph(Surv(t2,e2)~age+sex+sofa+charlson+lactate+creat+vlow,data=lm); b2<-summary(f2)$coefficients
cat(sprintf('  N=%d  VTI<17 HR=%.2f (%.2f-%.2f) p=%.3g\n', nrow(lm), exp(b2['vlow',1]),exp(b2['vlow',1]-1.96*b2['vlow',3]),exp(b2['vlow',1]+1.96*b2['vlow',3]),b2['vlow',5]))
cat('\n== 亚组 HR/CI ==\n')
for(s in list(c('通风','vent'),c('升压药',''))){ }
f3<-coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+vlow,data=dd[dd$vent=='t',]); b3<-summary(f3)$coefficients
cat(sprintf('  机械通气 (N=%d): HR=%.2f (%.2f-%.2f)\n', sum(dd$vent=='t'), exp(b3['vlow',1]),exp(b3['vlow',1]-1.96*b3['vlow',3]),exp(b3['vlow',1]+1.96*b3['vlow',3])))
