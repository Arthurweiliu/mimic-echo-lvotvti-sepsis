# 14_unit_of_analysis.R — unit of analysis: all eligible stays vs one record per patient
# Outputs: [A] all eligible ICU stays, [B] first ICU stay per patient (primary analysis),
#          [C] all stays with patient-clustered robust standard errors.
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
covs<-c('age','sex','sofa','charlson','lactate','creat','vti')
prep<-function(x){
  for(v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti')) x[[v]]<-num(x[[v]])
  x$sex<-factor(ifelse(x$sex %in% c('M','F'),x$sex,NA))
  x
}
suppressMessages(library(survival))
cov<-'age+sex+sofa+charlson+lactate+creat'

all<-prep(read.csv('h1_data_allstays.csv',stringsAsFactors=FALSE))     # every eligible stay
first<-prep(read.csv('h1_data.csv',stringsAsFactors=FALSE))           # one record per patient
cc<-function(x) complete.cases(x[,c('event90','fu_days',covs)]) & x$fu_days>0 & x$vti>0

a<-all[cc(all),]; a$vlow<-as.numeric(a$vti<17)
b<-first[cc(first),]; b$vlow<-as.numeric(b$vti<17)
cat('分析集(完整病例): 全部 stay N=',nrow(a),' 唯一患者=',length(unique(a$subject_id)),
    ' | 首次 ICU N=',nrow(b),' 唯一患者=',length(unique(b$subject_id)),'\n')
tb<-table(table(a$subject_id)); cat('每患者 stay 数分布:\n'); print(tb)
cat('>1 stay 的患者数=',sum(table(a$subject_id)>1),' 其贡献 stays=',sum(table(a$subject_id)[table(a$subject_id)>1]),'\n')

hr<-function(f,term='vlow'){ bb<-summary(f)$coefficients
  sprintf('%.2f (%.2f-%.2f) p=%.3g', exp(bb[term,1]), exp(bb[term,1]-1.96*bb[term,3]), exp(bb[term,1]+1.96*bb[term,3]), bb[term,5]) }
cat(sprintf('\n[A] 全部 stays            : N=%d  VTI<17 HR=%s\n', nrow(a), hr(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=a))))
cat(sprintf('[B] 仅首次 ICU(主分析)     : N=%d  VTI<17 HR=%s\n', nrow(b), hr(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=b))))
cat(sprintf('[C] 全部 stays + 患者聚类稳健: N=%d  VTI<17 HR=%s\n', nrow(a),
            hr(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow + cluster(subject_id)')),data=a))))
