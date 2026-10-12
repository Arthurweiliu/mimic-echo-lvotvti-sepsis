# 24_lvef_vs_vti.R — LVOT-VTI vs LVEF 头对头模型对比（+相关性）
setwd('F:/openclaw_workspace/sepsis_echo_lit')
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','charlson','lactate','creat','lvef','fu_days','event90','death28','vti')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA))
suppressMessages(library(survival))
dd<-d[complete.cases(d[,c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti','lvef')])&d$fu_days>0&d$vti>0,]
cat('N=',nrow(dd),' deaths=',sum(dd$event90),'\n')
cat('VTI-LVEF Spearman r=', round(cor(dd$vti, dd$lvef, method='spearman'),3),
    ' Pearson r=', round(cor(dd$vti, dd$lvef),3), '\n')
base<-'age+sex+sofa+charlson+lactate+creat'
c_ <- function(f) summary(coxph(as.formula(paste('Surv(fu_days,event90)~',f)),data=dd))$concordance[1]
cat('\n== C-index (90d) ==\n')
cat(' base                     :', round(c_(base),3), '\n')
cat(' base + LVEF              :', round(c_(paste(base,'+ lvef')),3), '\n')
cat(' base + VTI(low<17)       :', round(c_(paste(base,'+ I(vti<17)')),3), '\n')
cat(' base + VTI(continuous)   :', round(c_(paste(base,'+ vti')),3), '\n')
cat(' base + LVEF + VTI(low)   :', round(c_(paste(base,'+ lvef + I(vti<17)')),3), '\n')
auc<-function(f){ m<-glm(as.formula(paste('death28~',f)),data=dd,family=binomial); p<-predict(m,type='response'); y<-dd$death28; n1<-sum(y==1); n0<-sum(y==0); r<-rank(p); (sum(r[y==1])-n1*(n1+1)/2)/(n1*n0) }
cat('\n== 28d AUC ==\n')
cat(' base                     :', round(auc(base),3),'\n')
cat(' base + LVEF              :', round(auc(paste(base,'+ lvef')),3),'\n')
cat(' base + VTI(low<17)       :', round(auc(paste(base,'+ I(vti<17)')),3),'\n')
cat(' base + LVEF + VTI(low)   :', round(auc(paste(base,'+ lvef + I(vti<17)')),3),'\n')
