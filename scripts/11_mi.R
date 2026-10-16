# 11_mi.R — 多重插补（mice）敏感性：H1 主模型
# optional extra library path: set H1_LIB if your R packages live outside the default library
if (nzchar(Sys.getenv('H1_LIB'))) .libPaths(c(Sys.getenv('H1_LIB'), .libPaths()))
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','sapsii','charlson','lactate','creat','fu_days','event90','death28','vti')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA))
d<-d[!is.na(d$fu_days)&d$fu_days>0 & !is.na(d$event90),]
cat('可用行(结局齐全):',nrow(d),'\n')
suppressMessages({library(mice);library(survival)})
vars<-c('age','sex','sofa','charlson','lactate','creat','vti','fu_days','event90')
dd<-d[,vars]
cat('缺失: vti=',sum(is.na(dd$vti)),' creat=',sum(is.na(dd$creat)),' lactate=',sum(is.na(dd$lactate)),'\n')
set.seed(1); imp<-mice(dd, m=20, printFlag=FALSE)
poolcox<-function(expr_make){
  m<-imp$m; est<-se<-numeric(m)
  for(i in 1:m){ di<-complete(imp,i); di<-expr_make(di)
    f<-coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+vti_low, data=di)
    est[i]<-coef(f)['vti_low']; se[i]<-sqrt(vcov(f)['vti_low','vti_low']) }
  q<-mean(est); u<-mean(se^2); b<-var(est); T<-u+(1+1/m)*b; s<-sqrt(T)
  c(HR=exp(q), lo=exp(q-1.96*s), hi=exp(q+1.96*s), p=2*pnorm(-abs(q/s)))
}
cat('\n== 多重插补 (m=20)：VTI<17 vs >=17 (校正) ==\n')
r<-poolcox(function(di){di$vti_low<-as.numeric(di$vti<17); di}); cat(sprintf('  HR=%.2f (%.2f-%.2f) p=%.3g\n', r['HR'],r['lo'],r['hi'],r['p']))
# 连续
poolcont<-function(){
  m<-imp$m; est<-se<-numeric(m)
  for(i in 1:m){ di<-complete(imp,i); f<-coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+vti, data=di)
    est[i]<-coef(f)['vti']; se[i]<-sqrt(vcov(f)['vti','vti']) }
  q<-mean(est); u<-mean(se^2); b<-var(est); s<-sqrt(u+(1+1/m)*b)
  c(HR5=exp(q*5), lo=exp((q-1.96*s)*5), hi=exp((q+1.96*s)*5), p=2*pnorm(-abs(q/s)))
}
r2<-poolcont(); cat(sprintf('  每+5cm HR=%.2f (%.2f-%.2f) p=%.3g\n', r2['HR5'],r2['lo'],r2['hi'],r2['p']))
# complete-case comparison (same contrasts) computed here, not hard-coded
ccs<-d[complete.cases(d[,c('age','sex','sofa','charlson','lactate','creat','vti')]),]
fcc<-coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+I(vti<17), data=ccs); bc<-summary(fcc)$coefficients; ic<-grep('vti < 17',rownames(bc))
fcc2<-coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+vti, data=ccs); bc2<-summary(fcc2)$coefficients
cat(sprintf('\n(对照 完整病例 N=%d: <17 HR=%.2f %.2f-%.2f; 每+5cm HR=%.2f %.2f-%.2f)\n',
  nrow(ccs), exp(bc[ic,1]), exp(bc[ic,1]-1.96*bc[ic,3]), exp(bc[ic,1]+1.96*bc[ic,3]),
  exp(bc2['vti',1]*5), exp((bc2['vti',1]-1.96*bc2['vti',3])*5), exp((bc2['vti',1]+1.96*bc2['vti',3])*5)))
