# 18_km_tdhr.R — ① KM(含 number-at-risk，全部黑字) + 时间依赖 HR 曲线
# Style: all annotation text black; risk table has a labelled left gutter.
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA))
suppressMessages(library(survival))
dd<-d[complete.cases(d[,c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')])&d$fu_days>0&d$vti>0,]
dd$vcat2<-factor(ifelse(dd$vti<17,'VTI < 17 cm', ifelse(dd$vti<=20,'VTI 17-20 cm','VTI > 20 cm')),
                 levels=c('VTI < 17 cm','VTI 17-20 cm','VTI > 20 cm'))
dd$vlow<-as.numeric(dd$vti<17)
cols<-c('#c0392b','#e08a1e','#1f4e9c')          # print-safe version of red / orange / blue
lab <-c('< 17 cm','17-20 cm','> 20 cm')

# --- KM with number at risk ---
png('H1_KM_VTI.png', 3900, 3000, res=600)
par(mar=c(8.6,6.2,2.2,1.6))
fit<-survfit(Surv(fu_days,event90)~vcat2, data=dd)
sm<-summary(fit); ylo<-min(sm$lower)-0.03; if(ylo<0) ylo<-0
tt<-c(0,30,60,90)
plot(fit, col=cols, lwd=2.4, xlim=c(0,90), ylim=c(ylo,1.0), xaxt='n',
     xlab='', ylab='Survival probability', cex.lab=1.1, cex.axis=1.0)
axis(1, at=tt, labels=tt, line=0, cex.axis=1.0)
mtext('Days from ICU admission', side=1, line=1.9, cex=1.1)
legend('bottomleft', levels(dd$vcat2), col=cols, lwd=2.4, bty='n', cex=1.0)
mtext('Number at risk', side=1, line=3.7, adj=1, at=-1, cex=0.95, font=2)
for(i in 1:3){
  g<-levels(dd$vcat2)[i]
  nr<-sapply(tt, function(t) sum(dd$fu_days>=t & dd$vcat2==g))
  mtext(lab[i],            side=1, at=-1, adj=1, line=5.1+(i-1)*1.2, cex=0.9, col='black')
  mtext(as.character(nr),  side=1, at=tt, adj=0, line=5.1+(i-1)*1.2, cex=0.9, col='black')
}
dev.off(); cat('saved H1_KM_VTI.png (with number at risk)\n')

# --- time-dependent HR ---
f3<-coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+vlow+tt(vlow),
          tt=function(x,t,...) x*log(t+1), data=dd)
b<-coef(f3); V<-vcov(f3)
tv<-seq(1,90,length=200); lg<-log(tv+1)
lhr<-b['vlow']+b['tt(vlow)']*lg
se<-sqrt(V['vlow','vlow'] + lg^2*V['tt(vlow)','tt(vlow)'] + 2*lg*V['vlow','tt(vlow)'])
png('H1_tdHR.png', 3900, 3000, res=600)
par(mar=c(5.4,5.4,2.2,1.6))
yl<-range(exp(lhr-1.96*se), exp(lhr+1.96*se), 1)
plot(tv, exp(lhr), type='l', lwd=3, ylim=c(floor(yl[1]*10)/10, ceiling(yl[2]*10)/10),
     xlab='Days from ICU admission', ylab='Hazard ratio (VTI < 17 vs >= 17 cm)', cex.lab=1.1, cex.axis=1.0)
lines(tv, exp(lhr-1.96*se), lty=2, lwd=1.6); lines(tv, exp(lhr+1.96*se), lty=2, lwd=1.6)
abline(h=1, lty=3, col='grey40')
dev.off(); cat('saved H1_tdHR.png\n')
