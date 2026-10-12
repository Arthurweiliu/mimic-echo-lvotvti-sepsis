# 12_audit.R — 按协议反向核查研究过程
.libPaths(c('F:/openclaw_workspace/Rlibs', .libPaths()))
setwd('F:/openclaw_workspace/sepsis_echo_lit')
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','charlson','lactate','creat','fu_days','event90','death90','vti','norepi','mbp_min')) d[[v]]<-num(d[[v]])
cat('=== A. 入排核查 ===\n')
cat('N(队列)=',nrow(d),'\n')
cat('age 缺失=',sum(is.na(d$age)),' min=',min(d$age,na.rm=TRUE),' <18岁=',sum(d$age<18,na.rm=TRUE),'\n')
cat('sex 缺失/异常=',sum(!(d$sex %in% c('M','F'))),'\n')
cat('VTI 缺失(=NA或0)=',sum(is.na(d$vti)|d$vti==0),' 可得=',sum(!is.na(d$vti)&d$vti>0),'\n')
cat('VTI 极端: >40=',sum(d$vti>40,na.rm=TRUE),' <5=',sum(d$vti>0&d$vti<5,na.rm=TRUE),' >99? ',sum(d$vti>99,na.rm=TRUE),'\n')
cat('event90 缺失=',sum(is.na(d$event90)),' fu_days<=0=',sum(d$fu_days<=0,na.rm=TRUE),'\n')

cat('\n=== B. 时间原点与截尾 ===\n')
cat('fu_days 中位=',median(d$fu_days,na.rm=TRUE),' max=',max(d$fu_days,na.rm=TRUE),'(应<=90)\n')

cat('\n=== C. 主模型 PH 假设 (cox.zph) ===\n')
suppressMessages(library(survival))
dd<-d[complete.cases(d[,c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')])&d$fu_days>0&d$vti>0,]
m<-coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+I(vti<17),data=dd)
print(cox.zph(m))

cat('\n=== D. Table1 升压药异常复核 ===\n')
dd$vcat<-cut(dd$vti,c(-Inf,17,20,Inf),labels=c('<17','17-20','>20'))
cat('norepi 非缺失数=',sum(!is.na(dd$norepi)),' 其中>0=',sum(dd$norepi>0,na.rm=TRUE),'\n')
cat('把 NA 视为未用(0)时，各 VTI 组"用升压药%" :\n')
dd$vaso_any<-as.numeric(!is.na(dd$norepi)&dd$norepi>0)
for(g in levels(dd$vcat)) cat('  ',g,': ',round(mean(dd$vaso_any[dd$vcat==g])*100,1),'%\n')
cat('仅按非缺失计算(旧写法) 各 VTI 组% :\n')
for(g in levels(dd$vcat)){ x<-dd$norepi[dd$vcat==g]; cat('  ',g,': n_nonNA=',sum(!is.na(x)),' ',round(mean(x>0,na.rm=TRUE)*100,1),'%\n') }

cat('\n=== E. 时间验证分组 ===\n')
print(table(dd$yg))
cat('缺失 yg=',sum(is.na(dd$yg)|dd$yg==''),'\n')
