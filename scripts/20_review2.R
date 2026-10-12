# 20_review2.R — 二审核查计算
setwd('F:/openclaw_workspace/sepsis_echo_lit')
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','charlson','lactate','creat','lvef','fu_days','event90','vti','norepi')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA)); d$vent<-d$vent=='t'
library(survival)
cat('== 排除数算术 ==\n')
v<-!is.na(d$vti)&d$vti>0
cat('vti>0:',sum(v),'\n')
cc<-complete.cases(d[,c('age','sex','sofa','charlson','lactate','creat')])
cat('  vti>0 且协变量完整:',sum(v&cc),'\n')
cat('  vti>0 & 协变量完整 & fu_days>0:',sum(v&cc&d$fu_days>0),'\n')
cat('  无随访时间(fu<=0) 且 vti>0:',sum(v&d$fu_days<=0,na.rm=TRUE),'\n')
cat('  7,253-5,410 =',7253-5410,'\n')
cat('\n== 协变量缺失与无随访的构成 ==\n')
cat('  仅协变量缺失(vti>0, fu>0):',sum(v&!cc&d$fu_days>0),'\n')
cat('  仅无随访(vti>0, cc, fu<=0):',sum(v&cc&d$fu_days<=0),'\n')
cat('  两者皆(:',sum(v&!cc&d$fu_days<=0),'\n')

dd<-d[v&cc&d$fu_days>0,]
dd$vlow<-as.numeric(dd$vti<17)
cov<-'age+sex+sofa+charlson+lactate+creat'
cat('\n== LVEF 调整敏感性 ==\n')
base<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=dd); b<-summary(base)$coefficients
cat(sprintf('  主模型(<17): HR=%.2f (%.2f-%.2f)\n',exp(b['vlow',1]),exp(b['vlow',1]-1.96*b['vlow',3]),exp(b['vlow',1]+1.96*b['vlow',3])))
dl<-dd[!is.na(dd$lvef),]
f2<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ lvef + vlow')),data=dl); b2<-summary(f2)$coefficients
cat(sprintf('  加 LVEF(子集 N=%d): vlow HR=%.2f (%.2f-%.2f); LVEF HR/5%%=%.2f\n', nrow(dl),
  exp(b2['vlow',1]),exp(b2['vlow',1]-1.96*b2['vlow',3]),exp(b2['vlow',1]+1.96*b2['vlow',3]), exp(b2['lvef',1]*5)))
cat('\n== 亚组 HR/CI (VTI<17) ==\n')
for(s in c('vent','vaso')){
  sub<-if(s=='vent') dd[dd$vent,] else dd[!is.na(dd$norepi)&dd$norepi>0,]
  f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=sub); bb<-summary(f)$coefficients
  cat(sprintf('  %s (N=%d): HR=%.2f (%.2f-%.2f) p=%.3g\n', s, nrow(sub), exp(bb['vlow',1]),exp(bb['vlow',1]-1.96*bb['vlow',3]),exp(bb['vlow',1]+1.96*bb['vlow',3]),bb['vlow',5]))
}
