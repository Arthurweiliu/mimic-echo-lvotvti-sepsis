# 18_km_tdhr.R — ① KM(含 number-at-risk) + 时间依赖 HR 曲线
setwd('F:/openclaw_workspace/sepsis_echo_lit')
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA))
suppressMessages(library(survival))
dd<-d[complete.cases(d[,c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')])&d$fu_days>0&d$vti>0,]
dd$vcat2<-factor(ifelse(dd$vti<17,'VTI < 17 cm', ifelse(dd$vti<=20,'VTI 17-20 cm','VTI > 20 cm')),
                 levels=c('VTI < 17 cm','VTI 17-20 cm','VTI > 20 cm'))
dd$vlow<-as.numeric(dd$vti<17)
cols<-c('red','orange','blue')

# --- KM with number at risk ---
png('H1_KM_VTI.png', 3900, 3000, res=600)
par(mar=c(9,5.5,2,2))
fit<-survfit(Surv(fu_days,event90)~vcat2, data=dd)
plot(fit, col=cols, lwd=2, xlim=c(0,90), ylim=c(0.5,1.0), xaxt='n',
     xlab='Days from ICU admission', ylab='Survival probability')
legend('bottomleft', levels(dd$vcat2), col=cols, lwd=2, bty='n')
tt<-c(0,30,60,90)
mtext('Number at risk', side=1, line=3.4, adj=0, at=0, cex=0.85, font=2)
for(i in 1:3){
  g<-levels(dd$vcat2)[i]
  nr<-sapply(tt, function(t) sum(dd$fu_days>=t & dd$vcat2==g))
  mtext(as.character(nr), side=1, at=tt, line=4.4+(i-1)*1.1, col=cols[i], cex=0.8)
}
axis(1, at=tt, labels=tt, line=0)
dev.off(); cat('saved H1_KM_VTI.png (with number at risk)\n')

# --- time-dependent HR ---
f3<-coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+vlow+tt(vlow),
          tt=function(x,t,...) x*log(t+1), data=dd)
b<-coef(f3); V<-vcov(f3)
tv<-seq(1,90,length=200); lg<-log(tv+1)
lhr<-b['vlow']+b['tt(vlow)']*lg
se<-sqrt(V['vlow','vlow'] + lg^2*V['tt(vlow)','tt(vlow)'] + 2*lg*V['vlow','tt(vlow)'])
png('H1_tdHR.png', 3900, 3000, res=600)
par(mar=c(5.5,5.5,2,2))
plot(tv, exp(lhr), type='l', lwd=3, ylim=c(0.8, max(exp(lhr+1.96*se))+0.2),
     xlab='Days from ICU admission', ylab='Hazard ratio (VTI < 17 vs >= 17 cm)')
lines(tv, exp(lhr-1.96*se), lty=2, lwd=1.5); lines(tv, exp(lhr+1.96*se), lty=2, lwd=1.5)
abline(h=1, lty=3, col='grey40')
dev.off(); cat('saved H1_tdHR.png\n')
