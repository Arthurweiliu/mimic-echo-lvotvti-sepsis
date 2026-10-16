# 27_validity.R — prespecified "measurement-validity" analyses for LVOT-VTI
# The primary analysis is unchanged; this script adds:
#   (S9a) cohort restricted to patients in whom LVOT-VTI validly indexes forward stroke volume
#   (S9b) primary model additionally adjusted for AF/atrial arrhythmia, aortic valve disease, image quality
#   (S10) effect modification by AF, aortic stenosis >= moderate, image quality
#   (S11) quantitative bias analysis: E-value + regression-dilution scenarios for exposure measurement error
# Validity framework (see Methods): (A) velocity spectrum sampled at the subvalvular LVOT,
# (B) VTI x LVOT area = stroke volume (area assumed constant), (C) the recorded value represents
# a steady state (beat-to-beat averaging).
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x))
tr<-function(x) trimws(as.character(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti','rest_lvot','induc_lvot','lvot_pk_vel')) if(v %in% colnames(d)) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA))
suppressMessages(library(survival))

## ---- derived validity flags -------------------------------------------------
d$tamp_set   <- grepl('yes|early|collapse|tamponade', tr(d$tamponade), ignore.case=TRUE)
d$peff_mod   <- grepl('moderate|large|small-moderate', tr(d$pericardial_effusion), ignore.case=TRUE)
d$lvot_obstr <- (!is.na(d$rest_lvot) & d$rest_lvot>30) | (!is.na(d$induc_lvot) & d$induc_lvot>30) | (!is.na(d$lvot_pk_vel) & d$lvot_pk_vel>=2)
d$veg        <- grepl('mass|veg', tr(d$av_mass), ignore.case=TRUE)
d$as_mod     <- grepl('^mod|mod to severe|^severe|very severe|low flow', tr(d$aortic_stenosis), ignore.case=TRUE)
d$ar_mod     <- grepl('mod \\[2|mod-severe|severe \\[4', tr(d$aortic_regurg), ignore.case=TRUE)
d$tech       <- tr(d$technical_quality)
d$tech_poor  <- d$tech=='Poor'
d$tech_sub   <- d$tech %in% c('Suboptimal','Poor')
d$invalid_core <- d$tamp_set | d$peff_mod | d$lvot_obstr | d$veg
for(v in c('tamp_set','peff_mod','lvot_obstr','veg','as_mod','ar_mod','tech_poor','tech_sub','invalid_core')) d[[v]][is.na(d[[v]])]<-FALSE

## ---- atrial arrhythmia: ICU chart (Heart Rhythm) + ICD codes ----------------
af1<-read.csv('af_icu.csv',stringsAsFactors=FALSE); af2<-read.csv('af_icd.csv',stringsAsFactors=FALSE)
d$stay_id<-as.character(d$stay_id); af1$stay_id<-as.character(af1$stay_id); af2$stay_id<-as.character(af2$stay_id)
d<-merge(d, af1[,c('stay_id','af_icu')], by='stay_id', all.x=TRUE)
d<-merge(d, af2[,c('stay_id','af_icd')], by='stay_id', all.x=TRUE)
d$af<-as.numeric((!is.na(d$af_icu)&d$af_icu==1) | (!is.na(d$af_icd)&d$af_icd==1))

## ---- analysis set (identical to the primary analysis) ------------------------
ccs<-c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')
A<-d[complete.cases(d[,ccs]) & num(d$fu_days)>0 & num(d$vti)>0,]
A$vlow<-as.numeric(A$vti<17)
cov<-'age+sex+sofa+charlson+lactate+creat'
hr<-function(x,lab){ f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=x); b<-summary(f)$coefficients
  cat(sprintf('%-46s N=%5d  HR=%.2f (%.2f-%.2f) p=%.2g\n', lab, nrow(x), exp(b['vlow',1]), exp(b['vlow',1]-1.96*b['vlow',3]), exp(b['vlow',1]+1.96*b['vlow',3]), b['vlow',5])) }
cat('=== 队列与"测量有效性"限定集 (主分析 n =',nrow(A),') ===\n')
cat(sprintf('  心包填塞 %d | 中量以上积液 %d | LVOT梗阻 %d | 瓣周赘生物 %d | 技术质量Poor %d\n',
  sum(A$tamp_set), sum(A$peff_mod), sum(A$lvot_obstr), sum(A$veg), sum(A$tech_poor)))
cat(sprintf('  房颤/房扑(AF, ICU心律或ICD) n=%d (%.1f%%) | 主动脉瓣狭窄>=中度 %d | 反流>=中度 %d\n\n',
  sum(A$af), mean(A$af)*100, sum(A$as_mod), sum(A$ar_mod)))

cat('=== S9a. 限制在"测量有效性"人群 ===\n'); hr(A,'主分析(参考)')
hr(A[!A$invalid_core,], '排除核心无效集(压塞/积液/LVOT梗阻/赘生物)')
hr(A[!A$invalid_core & !A$tech_poor,], '  + 排除技术质量 Poor')
hr(A[!A$invalid_core & !A$tech_sub,], '  + 排除技术质量 Suboptimal/Poor')
hr(A[!(A$as_mod|A$ar_mod),], '排除主动脉瓣狭窄/反流 >=中度')
hr(A[!A$invalid_core & !(A$as_mod|A$ar_mod),], '核心无效集 + AS/AR>=中度 均排除')

cat('\n=== S9b. 主模型额外校正 ===\n')
b0<-summary(coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+vlow,data=A))$coefficients
b1<-summary(coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+vlow+af,data=A))$coefficients
b2<-summary(coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+vlow+af+as_mod+ar_mod,data=A))$coefficients
b3<-summary(coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+vlow+af+as_mod+ar_mod+tech_sub,data=A))$coefficients
for(nm in c('b0','b1','b2','b3')){ b<-get(nm)
  cat(sprintf('  %-28s HR(vlow)=%.2f (%.2f-%.2f) p=%.2g\n', nm, exp(b['vlow',1]), exp(b['vlow',1]-1.96*b['vlow',3]), exp(b['vlow',1]+1.96*b['vlow',3]), b['vlow',5])) }
cat(sprintf('  AF 本身: HR=%.2f (%.2f-%.2f) p=%.2g\n', exp(b3['af',1]), exp(b3['af',1]-1.96*b3['af',3]), exp(b3['af',1]+1.96*b3['af',3]), b3['af',5]))
cat(sprintf('  AS>=mod: HR=%.2f (%.2f-%.2f) | AR>=mod: HR=%.2f (%.2f-%.2f)\n',
  exp(b2['as_modTRUE',1]), exp(b2['as_modTRUE',1]-1.96*b2['as_modTRUE',3]), exp(b2['as_modTRUE',1]+1.96*b2['as_modTRUE',3]),
  exp(b2['ar_modTRUE',1]), exp(b2['ar_modTRUE',1]-1.96*b2['ar_modTRUE',3]), exp(b2['ar_modTRUE',1]+1.96*b2['ar_modTRUE',3])))

cat('\n=== S10. 效应修饰 (vlow x 修饰因素) ===\n')
for(m in c('af','as_mod','tech_sub')){
  f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow*',m)),data=A); b<-summary(f)$coefficients
  it<-grep(paste0('vlow:',m),rownames(b))
  cat(sprintf('  vlow x %-9s 交互 HR=%.2f (%.2f-%.2f) p=%.3f | 分层: %s\n', m,
    exp(b[it,1]), exp(b[it,1]-1.96*b[it,3]), exp(b[it,1]+1.96*b[it,3]), b[it,5],
    paste(sapply(c(0,1), function(g){ s<-A[A[[m]]==g,]; ff<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=s); bb<-summary(ff)$coefficients
      sprintf('%s: %.2f (%.2f-%.2f)', if(g==0)'no' else 'yes', exp(bb['vlow',1]), exp(bb['vlow',1]-1.96*bb['vlow',3]), exp(bb['vlow',1]+1.96*bb['vlow',3])) }), collapse=' | ')))
}

cat('\n=== S11. 定量偏倚分析 ===\n')
ev<-function(rr){ if(is.na(rr)) return(NA_real_); if(rr<1) rr<-1/rr; rr+sqrt(rr*(rr-1)) }
cat(sprintf('  E-value(主分析 HR=%.2f) = %.2f ; 校正后 HR=%.2f -> E-value = %.2f\n',
  exp(b0['vlow',1]), ev(exp(b0['vlow',1])), exp(b3['vlow',1]), ev(exp(b3['vlow',1]))))
bc<-summary(coxph(Surv(fu_days,event90)~age+sex+sofa+charlson+lactate+creat+vti,data=A))$coefficients
hr5<-exp(bc['vti',1]*5)
cat(sprintf('  连续暴露 per +5 cm HR=%.3f；按组内相关rho做回归稀释校正（情景）:\n', hr5))
for(rho in c(1,0.9,0.8,0.7)){
  cat(sprintf('    rho=%.1f -> 校正 HR=%.3f (%.3f-%.3f)\n', rho, exp(bc['vti',1]*5/rho),
    exp((bc['vti',1]-1.96*bc['vti',3])*5/rho), exp((bc['vti',1]+1.96*bc['vti',3])*5/rho))) }

## ---- write supplementary table ----------------------------------------------
f<-function(x) coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=x)
rows<-list()
add<-function(lab,x){ b<-summary(f(x))$coefficients
  data.frame(Analysis=lab, N=nrow(x), HR=round(exp(b['vlow',1]),2), CI_low=round(exp(b['vlow',1]-1.96*b['vlow',3]),2),
             CI_high=round(exp(b['vlow',1]+1.96*b['vlow',3]),2), p=format.pval(b['vlow',5],digits=2)) }
rows[[1]]<-add('Primary analysis (reference)',A)
rows[[2]]<-add('Validity-restricted: exclude tamponade / >=small-moderate effusion / LVOT obstruction / vegetation',A[!A$invalid_core,])
rows[[3]]<-add('  + exclude poor image quality',A[!A$invalid_core & !A$tech_poor,])
rows[[4]]<-add('  + exclude suboptimal or poor image quality',A[!A$invalid_core & !A$tech_sub,])
rows[[5]]<-add('Exclude aortic stenosis / regurgitation >= moderate',A[!(A$as_mod|A$ar_mod),])
out<-do.call(rbind,rows); write.csv(out,'TableS9_validity_restricted.csv',row.names=FALSE)

## S10 csv — effect modification
em<-list()
for(m in c('af','as_mod','tech_sub')){
  f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow*',m)),data=A); b<-summary(f)$coefficients
  it<-grep(paste0('vlow:',m),rownames(b))
  g<-lapply(c(0,1), function(k){ s<-A[A[[m]]==k,]; h<-summary(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow')),data=s))$coefficients
    c(n=nrow(s), HR=round(exp(h['vlow',1]),2), lo=round(exp(h['vlow',1]-1.96*h['vlow',3]),2), hi=round(exp(h['vlow',1]+1.96*h['vlow',3]),2)) })
  em[[length(em)+1]]<-data.frame(modifier=m, interaction_HR=round(exp(b[it,1]),2), int_lo=round(exp(b[it,1]-1.96*b[it,3]),2),
    int_hi=round(exp(b[it,1]+1.96*b[it,3]),2), int_p=format.pval(b[it,5],digits=2),
    stratum0_label=paste0(m,' absent'), stratum0_N=g[[1]]['n'], stratum0_HR=g[[1]]['HR'], stratum0_lo=g[[1]]['lo'], stratum0_hi=g[[1]]['hi'],
    stratum1_label=paste0(m,' present'), stratum1_N=g[[2]]['n'], stratum1_HR=g[[2]]['HR'], stratum1_lo=g[[2]]['lo'], stratum1_hi=g[[2]]['hi']) }
write.csv(do.call(rbind,em),'TableS10_effect_modification.csv',row.names=FALSE)

## S11 csv — bias analysis
bias<-data.frame(quantity=c('E-value, primary analysis (HR 1.38)','E-value, AF-adjusted (HR 1.34)',
  'per +5 cm, as observed (rho=1.0)','per +5 cm, dilution-corrected rho=0.9','per +5 cm, dilution-corrected rho=0.8','per +5 cm, dilution-corrected rho=0.7'),
  value=c(round(ev(exp(b0['vlow',1])),2), round(ev(exp(b3['vlow',1])),2), round(hr5,3),
          round(exp(bc['vti',1]*5/0.9),3), round(exp(bc['vti',1]*5/0.8),3), round(exp(bc['vti',1]*5/0.7),3)))
write.csv(bias,'TableS11_bias_analysis.csv',row.names=FALSE)
cat('saved TableS10_effect_modification.csv, TableS11_bias_analysis.csv\n')
cat('\nsaved TableS9_validity_restricted.csv\n'); print(out)
