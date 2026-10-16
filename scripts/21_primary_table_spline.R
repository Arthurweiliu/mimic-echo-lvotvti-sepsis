# 21_primary_table_spline.R — 主分析人群(5,410) Table1 + 限制性立方样条 HR 图
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','sapsii','charlson','lactate','creat','mbp_min','hr_mean','lvef','norepi',
           'bun','hgb','plt','wbc','bili','inr','vti','fu_days','event90')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA))
d$vent<-d$vent=='t'; d$rrt<-d$rrt=='t'; d$vaso<-as.numeric(!is.na(d$norepi)&d$norepi>0)
dd<-d[complete.cases(d[,c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')])&d$fu_days>0&d$vti>0,]
lev<-c('<17','17-20','>20'); dd$vcat<-factor(ifelse(dd$vti<17,'<17', ifelse(dd$vti<=20,'17-20','>20')), levels=lev)
med<-function(x) sprintf('%.1f (%.1f-%.1f)', median(x,na.rm=TRUE), quantile(x,.25,na.rm=TRUE), quantile(x,.75,na.rm=TRUE))
ppct<-function(x) sprintf('%d (%.0f%%)', sum(x,na.rm=TRUE), mean(x,na.rm=TRUE)*100)
pk<-function(v) tryCatch(format(kruskal.test(v~dd$vcat)$p.value,digits=2,scientific=TRUE),error=function(e)'NA')
pc<-function(v) tryCatch(format(suppressWarnings(chisq.test(table(dd$vcat,v)))$p.value,digits=2,scientific=TRUE),error=function(e)'NA')
rows<-list(); add<-function(lbl,f,pf) rows[[length(rows)+1]]<<-c(lbl, sapply(lev,function(g) f(dd[dd$vcat==g,])), pf)
add('N', function(x) as.character(nrow(x)),'')
add('Age, years', function(x) med(x$age), pk(dd$age))
add('Female, n (%)', function(x) sprintf('%d (%.0f%%)', sum(x$sex=='F'), mean(x$sex=='F')*100), pc(dd$sex=='F'))
add('SOFA (first day)', function(x) med(x$sofa), pk(dd$sofa))
add('SAPS-II', function(x) med(x$sapsii), pk(dd$sapsii))
add('Charlson index', function(x) med(x$charlson), pk(dd$charlson))
add('Lactate, mmol/L', function(x) med(x$lactate), pk(dd$lactate))
add('Creatinine, mg/dL', function(x) med(x$creat), pk(dd$creat))
add('Hemoglobin, g/dL', function(x) med(x$hgb), pk(dd$hgb))
add('Platelets, K/uL', function(x) med(x$plt), pk(dd$plt))
add('WBC, K/uL', function(x) med(x$wbc), pk(dd$wbc))
add('Bilirubin total, mg/dL', function(x) med(x$bili), pk(dd$bili))
add('INR', function(x) med(x$inr), pk(dd$inr))
add('Heart rate, mean, bpm', function(x) med(x$hr_mean), pk(dd$hr_mean))
add('LVEF, %', function(x) med(x$lvef), pk(dd$lvef))
add('Vasopressor (24 h), n (%)', function(x) sprintf('%d (%.0f%%)', sum(x$vaso), mean(x$vaso)*100), pc(dd$vaso))
add('Invasive MV (24 h), n (%)', function(x) ppct(x$vent), pc(dd$vent))
add('RRT (7 d), n (%)', function(x) ppct(x$rrt), pc(dd$rrt))
add('90-day mortality, n (%)', function(x) ppct(x$event90), pc(dd$event90))
tab<-do.call(rbind,rows); colnames(tab)<-c('Variable', lev, 'p')
write.csv(tab,'H1_Table1_primary.csv',row.names=FALSE)
cat('Primary-cohort Table1 N=',nrow(dd),' deaths=',sum(dd$event90),'\n'); print(tab)

# --- spline HR figure (reference VTI = 20 cm) ---
suppressMessages(library(survival))
fit<-coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+splines::ns(vti,3), data=dd)
vt<-seq(8,40,length=200)
nd<-data.frame(vti=vt, age=median(dd$age), sex=factor('F',levels=c('M','F')), sofa=median(dd$sofa),
               charlson=median(dd$charlson), lactate=median(dd$lactate), creat=median(dd$creat))
ref<-data.frame(vti=20, age=median(dd$age), sex=factor('F',levels=c('M','F')), sofa=median(dd$sofa),
                charlson=median(dd$charlson), lactate=median(dd$lactate), creat=median(dd$creat))
lp<-predict(fit, newdata=nd, type='lp'); lpr<-predict(fit, newdata=ref, type='lp')
X<-model.matrix(~splines::ns(vti,3), data=nd); Xr<-model.matrix(~splines::ns(vti,3), data=ref)
V<-vcov(fit)[grep('ns\\(vti',rownames(vcov(fit))), grep('ns\\(vti',rownames(vcov(fit)))]
dX<-sweep(X[, -1, drop=FALSE], 2, Xr[1,-1], '-')
lhr<-as.numeric(dX %*% coef(fit)[grep('ns\\(vti',names(coef(fit)))])
se<-sqrt(diag(dX %*% V %*% t(dX)))
xr <- c(8, 40)   # clinical/estimable range
lo <- exp(lhr-1.96*se); hi <- exp(lhr+1.96*se)
keep <- vt >= xr[1] & vt <= xr[2] & hi <= 4 & lo >= 0.3
ymax <- min(max(hi[keep]), 4); ymin <- max(min(lo[keep]), 0)
png('H1_spline_HR.png', 3900, 3000, res=600); par(mar=c(5.4,5.4,2.2,1.6))
plot(vt[keep], exp(lhr)[keep], type='l', lwd=3, ylim=c(ymin, ymax),
     xlab='LVOT-VTI (cm)', ylab='Adjusted HR (ref = 20 cm)', cex.lab=1.1, cex.axis=1.0)
lines(vt[keep], lo[keep], lty=2); lines(vt[keep], hi[keep], lty=2)
abline(h=1, lty=3, col='grey40'); abline(v=17, lty=4, col='black')
dev.off(); cat('saved H1_spline_HR.png\n')
