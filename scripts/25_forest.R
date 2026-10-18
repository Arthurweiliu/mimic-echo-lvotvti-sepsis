# 25_forest.R — Figure 3: adjusted hazard ratios for 90-day mortality (forest plot)
# Reconstructed as a dedicated, reproducible script (previously only an ad-hoc variant existed).
# Style: no embedded title, compact left margin, numeric HR (95% CI) column on the right.
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA))
suppressMessages(library(survival))
cov<-'age+sex+sofa+charlson+lactate+creat'
dd<-d[complete.cases(d[,c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')])&d$fu_days>0&d$vti>0,]
# same cut-point rule as 07_h1.R / 21_primary_table_spline.R: vti==17 and vti==20 fall in "17-20"
dd$vcat<-factor(ifelse(dd$vti<17,'<17', ifelse(dd$vti<=20,'17-20','>20')), levels=c('>20','<17','17-20'))

mc<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vcat')),data=dd); bc<-summary(mc)$coefficients
mn<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti')),data=dd);   bn<-summary(mn)$coefficients

mk <- function(lab, est, se){ hr<-exp(est); c(lab, hr, exp(est-1.96*se), exp(est+1.96*se)) }
items <- list(
  mk('LVOT-VTI < 17 vs > 20 cm',   bc['vcat<17',1],   bc['vcat<17',3]),
  mk('LVOT-VTI 17-20 vs > 20 cm',  bc['vcat17-20',1], bc['vcat17-20',3]),
  mk('LVOT-VTI per +5 cm',         bn['vti',1]*5,     bn['vti',3]*5))          # NOTE: *5 (per 5 cm)

k<-length(items); ys<-k:1
xlo<-0.85; xhi<-1.60; p1<-4.00; p2<-7.45
xp<-function(h) p1 + (h-xlo)/(xhi-xlo)*(p2-p1)

png('H1_forest.png', 3900, 2400, res=600)
par(mar=c(5.2,0.6,2.4,0.6))
plot(NA, xlim=c(0,10), ylim=c(0.3,k+0.7), axes=FALSE, xlab='', ylab='')
rect(p1, 0.3, p2, k+0.7, border='grey70', lwd=1)
abline(v=xp(1), lty=2, col='grey40')
for(i in 1:k){
  hr<-as.numeric(items[[i]][2]); lo<-as.numeric(items[[i]][3]); hi<-as.numeric(items[[i]][4])
  segments(xp(lo), ys[i], xp(hi), ys[i], lwd=2.4)
  points(xp(hr), ys[i], pch=18, cex=1.6)
  text(p1-0.30, ys[i], items[[i]][1], adj=1, cex=1.0)
  text(p2+0.30, ys[i], sprintf('%.2f (%.2f–%.2f)', hr, lo, hi), adj=0, cex=1.0)
}
ticks<-c(0.9,1.0,1.1,1.2,1.3,1.4,1.5)
axis(1, at=xp(ticks), labels=sprintf('%.1f',ticks), cex.axis=0.92, gap.axis=0)
mtext('Adjusted hazard ratio (90-day mortality)', side=1, line=2.6, cex=1.1)
text(p2+0.30, k+0.55, 'HR (95% CI)', adj=0, cex=1.0, font=2)
dev.off(); cat('saved H1_forest.png\n')
