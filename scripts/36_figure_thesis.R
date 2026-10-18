# 36_figure_thesis.R — core figure for the "per-beat flow vs cardiac output" thesis
{ .a <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
d <- read.csv('h1_data.csv', stringsAsFactors=FALSE)
for (v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti','vti_diam','lvef')) d[[v]]<-num(d[[v]])
d$sex <- factor(ifelse(d$sex %in% c('M','F'), d$sex, NA)); d$stay_id <- as.character(d$stay_id)
suppressMessages(library(survival))
cov <- 'age+sex+sofa+charlson+lactate+creat'
ad0 <- read.csv('analysis_data.csv', stringsAsFactors=FALSE, colClasses=c(stay_id='character'))
ad0 <- ad0[!duplicated(ad0$stay_id), c('stay_id','measurement_id','body_surface_area')]
d <- merge(d, ad0, by='stay_id', all.x=TRUE); d$bsa<-num(d$body_surface_area); d$measurement_id<-as.character(d$measurement_id)
ex <- read.csv('echo_extra.csv', stringsAsFactors=FALSE); ex$measurement_id<-as.character(ex$measurement_id)
hr <- ex[ex$measurement=='resting_hr', c('measurement_id','result')]; hr<-hr[!duplicated(hr$measurement_id),]; names(hr)[2]<-'rhr'
d <- merge(d, hr, by='measurement_id', all.x=TRUE); d$rhr<-num(d$rhr)
ccs <- c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')
A <- d[complete.cases(d[,ccs]) & num(d$fu_days)>0 & num(d$vti)>0, ]
A$sv <- A$vti*pi*(A$vti_diam/2)^2
A$hr <- ifelse(!is.na(A$rhr) & A$rhr>20 & A$rhr<220, A$rhr, NA)
A$co <- A$sv*A$hr/1000
A$ci <- ifelse(!is.na(A$bsa) & A$bsa>0.5, A$co/A$bsa, NA)
H <- A[!is.na(A$sv)&A$sv>5&A$sv<250 & !is.na(A$co)&A$co>1&A$co<15 & !is.na(A$ci)&A$ci>0.5&A$ci<8, ]
H$vlow <- as.numeric(H$vti<17); H$tachy <- as.numeric(H$hr>100)

expos <- list(c('vti',5,'LVOT-VTI (+5 cm)'), c('sv',10,'Stroke volume (+10 mL)'),
  c('lvef',10,'LVEF (+10 %)'), c('co',1,'Cardiac output (+1 L/min)'),
  c('ci',1,'Cardiac index (+1 L/min/m²)'), c('hr',10,'Heart rate (+10 bpm)'))
res <- data.frame()
for (e in expos){ var<-e[1]; k<-as.numeric(e[2]); lab<-e[3]
  f <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+',var)), data=H); b<-summary(f)$coefficients
  res <- rbind(res, data.frame(lab=lab, est=exp(b[var,1]*k), lo=exp((b[var,1]-1.96*b[var,3])*k),
    hi=exp((b[var,1]+1.96*b[var,3])*k), p=b[var,5]))
}
res$lab <- factor(res$lab, levels=rev(res$lab))
write.csv(res,'Item4_thesis_forest.csv', row.names=FALSE)
print(res)

png('Item4_thesis_figure.png', width=3900, height=2900, res=600)
layout(matrix(c(1,1,2,2,3,3), ncol=2, byrow=TRUE), heights=c(1.25,1,1))
# Panel A forest
par(mar=c(4,13,3,2))
xr <- range(c(res$lo,res$hi)); xr <- c(min(xr,0.85), max(xr,1.15))
plot(NA, xlim=log(xr), ylim=c(0.5,nrow(res)+0.5), yaxt='n', xaxt='n', xlab='Adjusted hazard ratio (95% CI)', ylab='',
     main='A. What carries prognostic information in sepsis?', cex.main=1.15, col.main='black')
abline(v=0, lty=2, col='grey60')
cols <- ifelse(grepl('Heart rate', res$lab), '#b2182b', ifelse(grepl('output|index', res$lab), 'grey45', '#2166ac'))
for (i in 1:nrow(res)){ y<-i; segments(log(res$lo[i]),y,log(res$hi[i]),y,col=cols[i],lwd=2.6)
  points(log(res$est[i]),y,pch=19,col=cols[i],cex=1.3) }
axis(1, at=log(c(0.85,0.9,0.95,1.0,1.05,1.1)), labels=c('0.85','0.90','0.95','1.00','1.05','1.10'))
axis(2, at=1:nrow(res), labels=res$lab, las=1, cex.axis=0.95)
# Panel B interaction by tachycardia
par(mar=c(4,9.5,3,2))
strat <- data.frame()
for (v in c('vlow','co')){ for (g in 0:1){ s<-H[H$tachy==g,]
  f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+',v)),data=s); b<-summary(f)$coefficients
  strat<-rbind(strat, data.frame(expo=ifelse(v=='vlow','VTI <17 (per-beat)','CO (per-min)'), grp=ifelse(g==0,'HR <=100','HR >100'),
    est=ifelse(v=='vlow',exp(b[v,1]),exp(b[v,1]*1)), lo=ifelse(v=='vlow',exp(b[v,1]-1.96*b[v,3]),exp((b[v,1]-1.96*b[v,3]))),
    hi=ifelse(v=='vlow',exp(b[v,1]+1.96*b[v,3]),exp((b[v,1]+1.96*b[v,3]))))) } }
print(strat)
plot(NA, xlim=c(0.7,1.6), ylim=c(0.5,4.5), yaxt='n', xaxt='n', ylab='', xlab='Adjusted HR', main='B. Why: HR is built into CO', cex.main=1.15)
abline(v=1, lty=2, col='grey60')
axis(1, at=c(0.7,0.9,1.0,1.2,1.4)); 
ys <- c(4,3,1.5,0.5)
labs <- c('VTI<17, HR<=100','VTI<17, HR>100','CO, HR<=100','CO, HR>100')
for (i in 1:4){ s<-strat[i,]; segments(s$lo,ys[i],s$hi,ys[i],col=c('#2166ac','#2166ac','grey45','grey45')[i],lwd=2.6)
  points(s$est,ys[i],pch=19,col=c('#2166ac','#2166ac','grey45','grey45')[i],cex=1.3) }
axis(2, at=ys, labels=labs, las=1, cex.axis=0.85)
# Panel C schematic
par(mar=c(1,2,1,2)); plot.new()
text(0.5,0.86,'C. The cancellation', cex=1.15, font=2)
text(0.5,0.62,'SV (protective)  x  HR (adverse)  =  CO (no association)', cex=1.0)
text(0.5,0.42,'HR 1.09 / 10 bpm   |   SV 0.94 / 10 mL   ->   CO 0.99 / 1 L/min (NS)', cex=0.85, col='grey25')
text(0.5,0.22,'A per-beat measure keeps the signal; a per-minute measure cancels it.', cex=0.9, font=3, col='grey20')
dev.off()
cat('\nsaved Item4_thesis_figure.png\n'); print(strat)
