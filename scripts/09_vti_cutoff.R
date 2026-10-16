# 09_vti_cutoff.R — VTI 切点依据：扫描 + 样条 LRT + 图
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA))
suppressMessages(library(survival))
cov<-'age+sex+sofa+charlson+lactate+creat'
cc<-complete.cases(d[,c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')])&d$fu_days>0&d$vti>0
dd<-d[cc,]
cat('N=',nrow(dd),' 事件=',sum(dd$event90==1),'\n\n')
cps<-13:26; ps<-numeric(length(cps)); best<-list(cp=NA,p=Inf)
for(i in seq_along(cps)){ cp<-cps[i]; dd$low<-as.numeric(dd$vti<cp)
  m<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ low')),data=dd)
  ps[i]<-summary(m)$coefficients['low',5]; if(ps[i]<best$p) best<-list(cp=cp,p=ps[i]) }
cat('== 切点扫描 (Cox LRT) ==\n'); for(i in seq_along(cps)) cat(sprintf('  <%-2d cm: p=%.2e%s\n', cps[i], ps[i], if(cps[i]==17)'   <- 文献 17' else ''))
cat('  数据最优 ≈', best$cp, 'cm (p=', format(best$p,digits=3), ')\n')
m0<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov)),data=dd)
m4<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ splines::ns(vti,3)')),data=dd)
a<-anova(m0,m4,test='LRT'); pspline<-a[2, ncol(a)]
cat('  样条(3df) LRT p =', format(pspline,digits=3), '\n')
cat('  参考 Rolston 2026 框架: <17 低 / 17-20 正常 / >20 高动力\n')
png('H1_cutoff_scan.png',3900,2700,res=600)
par(mar=c(5.4,5.4,2.6,1.6))
ymax<-max(-log10(ps))
plot(cps, -log10(ps), type='b', pch=19, xlab='VTI cutpoint (cm)', ylab='-log10(LRT p)',
     ylim=c(0, ymax*1.18), cex.lab=1.1, cex.axis=1.0)
abline(v=17, col='red', lty=2, lwd=2)
abline(v=best$cp, col='blue', lty=3, lwd=2)
text(best$cp+0.15, ymax*1.14, paste0('data-optimal ', best$cp,' cm'), col='black', pos=4, cex=0.95)
text(17.15, ymax*1.01, 'literature 17 cm', col='black', pos=4, cex=0.95)
dev.off(); cat('saved H1_cutoff_scan.png\n')
