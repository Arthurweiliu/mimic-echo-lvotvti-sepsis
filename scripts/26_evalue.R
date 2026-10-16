# 26_evalue.R — Table S5: E-value for unmeasured confounding (primary contrast)
# E-value for a hazard-ratio estimate treated as a risk ratio (VanderWeele & Ding 2017):
#   E = RR + sqrt(RR*(RR-1))     [RR >= 1]
# Computed for the primary estimate and for its lower 95% confidence limit.
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA))
suppressMessages(library(survival))
dd<-d[complete.cases(d[,c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')])&d$fu_days>0&d$vti>0,]
# same cut-point rule as 07_h1.R/21_primary_table_spline.R: vti==17 and vti==20 fall in "17-20"
dd$vcat<-factor(ifelse(dd$vti<17,'<17', ifelse(dd$vti<=20,'17-20','>20')), levels=c('>20','<17','17-20'))
m<-coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+vcat, data=dd)
b<-summary(m)$coefficients; ci<-summary(m)$conf.int
hr<-ci['vcat<17',1]; lo<-ci['vcat<17',3]
evalue<-function(rr){ if(is.na(rr)) return(NA_real_); if(rr<1) rr<-1/rr; rr+sqrt(rr*(rr-1)) }
ev<-evalue(hr); evlo<-evalue(lo)
cat(sprintf('Primary contrast (VTI < 17 vs > 20 cm): HR=%.3f (%.3f-%.3f)\n', hr, lo, ci['vcat<17',4]))
cat(sprintf('E-value, point estimate = %.2f ; lower 95%% CI = %.2f\n', ev, evlo))
out<-data.frame(contrast='LVOT-VTI < 17 vs > 20 cm', HR=round(hr,3), CI_low=round(lo,3), CI_high=round(ci['vcat<17',4],3),
                E_value=round(ev,2), E_value_CI=round(evlo,2), N=nrow(dd))
write.csv(out,'TableS5_evalue.csv',row.names=FALSE)
cat('saved TableS5_evalue.csv\n')
