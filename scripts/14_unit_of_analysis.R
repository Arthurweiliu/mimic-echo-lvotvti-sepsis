# 14_unit_of_analysis.R — 分析单元核查：多 ICU stay/患者的处理
setwd('F:/openclaw_workspace/sepsis_echo_lit')
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
cc<-read.csv('cohort_covariates.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA))
m<-merge(d, cc[,c('stay_id','icu_intime','first_icu_stay')], by='stay_id', all.x=TRUE)
suppressMessages(library(survival))
dd<-m[complete.cases(m[,c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')])&m$fu_days>0&m$vti>0,]
cat('分析集(完整病例) N stays=',nrow(dd),' 唯一患者=',length(unique(dd$subject_id)),'\n')
tab<-table(table(dd$subject_id)); cat('每患者 stay 数分布:\n'); print(tab)
cat('>1 stay 的患者数=',sum(table(dd$subject_id)>1),' 其贡献 stays=',sum(table(dd$subject_id)[table(dd$subject_id)>1]),'\n')
dd$vlow<-as.numeric(dd$vti<17)
cov<-'age+sex+sofa+charlson+lactate+creat'
f_all<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=dd)
b<-summary(f_all)$coefficients
cat(sprintf('\n[A] 全部 stays:  VTI<17 HR=%.2f (%.2f-%.2f) p=%.3g\n', exp(b['vlow',1]),exp(b['vlow',1]-1.96*b['vlow',3]),exp(b['vlow',1]+1.96*b['vlow',3]),b['vlow',5]))
# [B] 每患者仅首次 ICU（最早 intime）
m2<-m[complete.cases(m[,c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti','icu_intime')])&m$fu_days>0&m$vti>0,]
m2<-m2[order(m2$subject_id, m2$icu_intime),]
first<-m2[!duplicated(m2$subject_id),]
first$vlow<-as.numeric(first$vti<17)
f_f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=first); bf<-summary(f_f)$coefficients
cat(sprintf('[B] 仅首次 ICU: N=%d (患者=%d) VTI<17 HR=%.2f (%.2f-%.2f) p=%.3g\n', nrow(first), length(unique(first$subject_id)), exp(bf['vlow',1]),exp(bf['vlow',1]-1.96*bf['vlow',3]),exp(bf['vlow',1]+1.96*bf['vlow',3]),bf['vlow',5]))
# [C] 全部 stays + 患者聚类稳健方差
f_c<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow + cluster(subject_id)')),data=dd); bc<-summary(f_c)$coefficients
cat(sprintf('[C] 全部+聚类稳健: VTI<17 HR=%.2f (%.2f-%.2f) p=%.3g\n', exp(bc['vlow',1]),exp(bc['vlow',1]-1.96*bc['vlow',3]),exp(bc['vlow',1]+1.96*bc['vlow',3]),bc['vlow',5]))
