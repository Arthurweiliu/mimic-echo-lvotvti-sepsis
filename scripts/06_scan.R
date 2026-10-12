# 06_scan.R — 候选假设快速信号扫描（base R + survival）
setwd('F:/openclaw_workspace/sepsis_echo_lit')
num <- function(x) suppressWarnings(as.numeric(x))
a <- read.csv('analysis_data.csv', stringsAsFactors=FALSE)
c <- read.csv('congestion_analysis.csv', stringsAsFactors=FALSE)
for (v in c('lvot_vti','lvot_diam','mean_E_to_eprime','e_e_prime','tr_velocity','la_dimen','lvef')) a[[v]] <- num(a[[v]])
m <- merge(c, a[,c('stay_id','lvot_vti','lvot_diam','mean_E_to_eprime','e_e_prime','rv_function')], by='stay_id', all.x=TRUE)
for (v in c('age','sofa','sapsii','charlson','lactate','creat','fu_days','event90','death28','death90','aki_any','aki3','crrt','score3','RVabn','highTR','LAabn')) m[[v]] <- num(m[[v]])
m$sex <- factor(m$sex)
m$vent <- factor(m$vent24)
suppressMessages(library(survival))
cov <- 'age+sex+sofa+charlson+lactate+creat'
HR <- function(fit, term){ b <- summary(fit)$coefficients
  if(!(term %in% rownames(b))) return('term not found')
  sprintf('HR=%.2f (%.2f-%.2f) p=%.3g', exp(b[term,1]), exp(b[term,1]-1.96*b[term,3]), exp(b[term,1]+1.96*b[term,3]), b[term,5]) }

cat('###### H1: LVOT VTI（低每搏量）→ 90d 死亡 ######\n')
dd <- m[complete.cases(m[,c('event90','fu_days','lvot_vti','age','sex','sofa','charlson','lactate','creat')]) & m$fu_days>0, ]
cat('N=',nrow(dd),' 事件=',sum(dd$event90==1),'\n')
dd$lowVTI <- as.numeric(dd$lvot_vti < 18)
cat('低VTI(<18cm)比例=',round(mean(dd$lowVTI)*100,1),'%\n')
cat('  低VTI vs 非低 90d死亡率:', round(tapply(dd$event90, dd$lowVTI, mean)*100,1), '\n')
f <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ lowVTI')), data=dd); cat('  lowVTI:', HR(f,'lowVTI'), '\n')
f2 <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ lvot_vti')), data=dd); cat('  per +5cm VTI:', HR(f2,'lvot_vti'), '\n')

cat('\n###### H2: LVEF 分型（含高动力>70%）→ 90d 死亡 ######\n')
d2 <- m[complete.cases(m[,c('event90','fu_days','lvef','age','sex','sofa','charlson','lactate','creat')]) & m$fu_days>0, ]
d2$lvefcat <- cut(d2$lvef, c(-Inf,40,70,Inf), labels=c('lt40','40to70','gt70'))
cat('N=',nrow(d2),'\n'); cat('  90d死亡率按LVEF:', round(tapply(d2$event90, d2$lvefcat, mean)*100,2), '\n')
f <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ lvefcat')), data=d2); cat('  ', HR(f,'lvefcatlt40'), '\n  ', HR(f,'lvefcatgt70'), '\n')

cat('\n###### H3: 高 E/e′（充盈压）→ 90d 死亡 ######\n')
d3 <- m[complete.cases(m[,c('event90','fu_days','mean_E_to_eprime','age','sex','sofa','charlson','lactate','creat')]) & m$fu_days>0, ]
cat('N=',nrow(d3),' 事件=',sum(d3$event90==1),'\n')
d3$hiEe <- as.numeric(d3$mean_E_to_eprime >= 15)
cat('  高E/e′(≥15)比例=',round(mean(d3$hiEe)*100,1),'%  死亡率:', round(tapply(d3$event90, d3$hiEe, mean)*100,1), '\n')
f <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ hiEe')), data=d3); cat('  hiEe:', HR(f,'hiEe'), '\n')
f2 <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ mean_E_to_eprime')), data=d3); cat('  per +5 E/e′:', HR(f2,'mean_E_to_eprime'), '\n')

cat('\n###### H4: 右室异常 × 机械通气 交互 → 90d 死亡 ######\n')
d4 <- m[complete.cases(m[,c('event90','fu_days','RVabn','vent','age','sex','sofa','charlson','lactate','creat')]) & m$fu_days>0, ]
cat('N=',nrow(d4),'\n')
f <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ RVabn*vent')), data=d4)
b<-summary(f)$coefficients
cat('  RVabn:', HR(f,'RVabn'), '\n'); cat('  ventt:', HR(f,'ventt'), '\n'); cat('  交互 RVabn:ventt:', HR(f,'RVabn:ventt'), '\n')

cat('\n###### H5: 淤血评分 × 血管活性药 交互 → 90d 死亡 ######\n')
d5 <- m[complete.cases(m[,c('event90','fu_days','score3','norepi','age','sex','sofa','charlson','lactate','creat')]) & m$fu_days>0, ]
d5$hiVaso <- as.numeric(d5$norepi > 0)
f <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ score3*hiVaso')), data=d5)
cat('N=',nrow(d5),'\n'); cat('  score3:', HR(f,'score3'), '\n'); cat('  交互 score3:hiVaso:', HR(f,'score3:hiVaso'), '\n')
