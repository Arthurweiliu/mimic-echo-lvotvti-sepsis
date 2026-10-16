# 10_h1_combo.R — 第1步：VTI 低 × 淤血高 联合分层；+ 缺失/敏感性
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
h<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
c<-read.csv('congestion_analysis.csv',stringsAsFactors=FALSE)
m<-merge(h, c[,c('stay_id','score3','highTR','LAabn','RVabn','highIVC')], by='stay_id', all.x=TRUE)
for(v in c('age','sofa','charlson','lactate','creat','fu_days','event90','death28','vti','score3','highTR','LAabn','RVabn','highIVC')) m[[v]]<-num(m[[v]])
m$sex<-factor(ifelse(m$sex %in% c('M','F'),m$sex,NA))
suppressMessages(library(survival))
HRr<-function(fit,term){ b<-summary(fit)$coefficients
  if(!(term %in% rownames(b))) return('NA')
  sprintf('%.2f (%.2f-%.2f) p=%.2g', exp(b[term,1]), exp(b[term,1]-1.96*b[term,3]), exp(b[term,1]+1.96*b[term,3]), b[term,5]) }
HRr5<-function(fit,term){ b<-summary(fit)$coefficients      # per +5 cm (coef is per +1 cm)
  if(!(term %in% rownames(b))) return('NA')
  sprintf('%.2f (%.2f-%.2f) p=%.2g', exp(b[term,1]*5), exp((b[term,1]-1.96*b[term,3])*5), exp((b[term,1]+1.96*b[term,3])*5), b[term,5]) }
cov<-'age+sex+sofa+charlson+lactate+creat'

cc<-complete.cases(m[,c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti','score3')])&m$fu_days>0&m$vti>0
d<-m[cc,]; d$congHigh<-as.numeric(d$score3>=2); d$vtiLow<-as.numeric(d$vti<17)
d$grp<-factor(paste0(ifelse(d$vtiLow==1,'VTIlo','VTIhi'),'_',ifelse(d$congHigh==1,'CongHi','CongLo')),
              levels=c('VTIhi_CongLo','VTIhi_CongHi','VTIlo_CongLo','VTIlo_CongHi'))
cat('== 分组 N 与 90d 死亡率 ==\n')
for(g in levels(d$grp)){ x<-d$event90[d$grp==g]; cat(sprintf('  %-14s N=%4d  90d=%.1f%%\n', g, length(x), mean(x)*100)) }
f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ grp')),data=d)
cat('\n== 校正 HR（参照 VTIhi_CongLo）==\n')
for(g in levels(d$grp)[-1]) cat('  ',g,':',HRr(f,paste0('grp',g)),'\n')
f2<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vtiLow*congHigh')),data=d)
cat('\n== 交互 vtiLow:congHigh ==\n  ',HRr(f2,'vtiLow:congHigh'),'\n')
cat('  vtiLow:',HRr(f2,'vtiLow'),'| congHigh:',HRr(f2,'congHigh'),'\n')

cat('\n== 缺失分析：VTI 可得 vs 不可得 ==\n')
m$vtiAvail<-as.numeric(!is.na(m$vti)&m$vti>0)
for(v in c('age','sofa','sapsii','charlson','lactate','creat','death90')) {
  a<-mean(m[[v]][m$vtiAvail==1],na.rm=TRUE); b<-mean(m[[v]][m$vtiAvail==0],na.rm=TRUE)
  cat(sprintf('  %-8s VTI有=%.1f  VTI无=%.1f\n', v, a, b))
}
cat('\n== 敏感性 ==\n')
s1<-d[!(d$fu_days<=1 & d$event90==1),]; f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti')),data=s1)
cat('  排除<24h早死: per+5cm',HRr5(f,'vti'),' (N=',nrow(s1),')\n')
s2<-d[d$vti<=40,]; f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti')),data=s2)
cat('  排除VTI>40: per+5cm',HRr5(f,'vti'),' (N=',nrow(s2),')\n')
s3<-d[d$vti>=5,]; f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti')),data=s3)
cat('  排除VTI<5 : per+5cm',HRr5(f,'vti'),' (N=',nrow(s3),')\n')
