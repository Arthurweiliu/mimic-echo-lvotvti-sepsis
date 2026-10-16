# 38_fdr_panel.R — exact BH-FDR for the effect-modification family (Set A)
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=','',.a[1]))) else getwd()) }
num<-function(x) suppressWarnings(as.numeric(x)); tr<-function(x) trimws(as.character(x))
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
for(v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti','vti_diam','hr_mean')) d[[v]]<-num(d[[v]])
d$sex<-factor(ifelse(d$sex %in% c('M','F'),d$sex,NA)); d$stay_id<-as.character(d$stay_id)
suppressMessages(library(survival)); cov<-'age+sex+sofa+charlson+lactate+creat'
ad<-read.csv('analysis_data.csv',stringsAsFactors=FALSE,colClasses=c(stay_id='character'))
ad<-ad[!duplicated(ad$stay_id),c('stay_id','measurement_id')]
d<-merge(d,ad,by='stay_id',all.x=TRUE); d$measurement_id<-as.character(d$measurement_id)
ex<-read.csv('echo_extra.csv',stringsAsFactors=FALSE); ex$measurement_id<-as.character(ex$measurement_id)
pick<-function(m,nm){s<-ex[ex$measurement==m,c('measurement_id','result')];s<-s[!duplicated(s$measurement_id),];names(s)[2]<-nm;s}
d<-merge(d,pick('resting_hr','rhr'),by='measurement_id',all.x=TRUE)
d<-merge(d,pick('mitral_regurg','mr'),by='measurement_id',all.x=TRUE)
d<-merge(d,pick('tricuspid_regurg','trr'),by='measurement_id',all.x=TRUE)
af1<-read.csv('af_icu.csv',stringsAsFactors=FALSE);af2<-read.csv('af_icd.csv',stringsAsFactors=FALSE)
af1$stay_id<-as.character(af1$stay_id);af2$stay_id<-as.character(af2$stay_id)
d<-merge(d,af1[,c('stay_id','af_icu')],by='stay_id',all.x=TRUE);d<-merge(d,af2[,c('stay_id','af_icd')],by='stay_id',all.x=TRUE)
d$af<-as.numeric((!is.na(d$af_icu)&d$af_icu==1)|(!is.na(d$af_icd)&d$af_icd==1)); d$rhr<-num(d$rhr)
d$as_mod<-grepl('^mod|mod to severe|^severe|very severe|low flow',tr(d$aortic_stenosis),ignore.case=TRUE)
d$ar_mod<-grepl('mod \\[2|mod-severe|severe \\[4',tr(d$aortic_regurg),ignore.case=TRUE)
d$mr_mod<-grepl('mod',tr(d$mr),ignore.case=TRUE)&!grepl('mild',tr(d$mr),ignore.case=TRUE)
d$tr_mod<-grepl('mod|severe',tr(d$trr),ignore.case=TRUE)&!grepl('mild|trace',tr(d$trr),ignore.case=TRUE)
d$tech_sub<-tr(d$technical_quality)%in%c('Suboptimal','Poor')
d$tachyE<-(!is.na(d$rhr)&d$rhr>100); d$tachyI<-(!is.na(d$hr_mean)&d$hr_mean>100)
for(v in c('as_mod','ar_mod','mr_mod','tr_mod','tech_sub','tachyE','tachyI')) d[[v]][is.na(d[[v]])]<-FALSE
cc<-c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')
A<-d[complete.cases(d[,cc])&num(d$fu_days)>0&num(d$vti)>0,]; A$vlow<-as.numeric(A$vti<17)
cat('Set A N =',nrow(A),'\n')
ip<-function(m){f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vlow*',m)),data=A);b<-summary(f)$coefficients
  it<-grep(paste0('^vlow:',m,'TRUE$|^vlow:',m,'$'),rownames(b))[1]
  c(p=b[it,5], hr=exp(b[it,1]), lo=exp(b[it,1]-1.96*b[it,3]), hi=exp(b[it,1]+1.96*b[it,3]))}
show<-function(mods,lab){r<-do.call(rbind,lapply(mods,ip)); r<-data.frame(modifier=mods,r)
  r$FDR<-round(p.adjust(r$p,method='BH'),3); r$p<-signif(r$p,2)
  cat('\n---',lab,'---\n'); print(r,row.names=FALSE,digits=3)}
show(c('af','tachyE','as_mod','ar_mod','mr_mod','tr_mod','tech_sub'),'7-test family, tachycardia = ECHO HR')
show(c('af','tachyI','as_mod','ar_mod','mr_mod','tr_mod','tech_sub'),'7-test family, tachycardia = ICU HR')
show(c('af','as_mod','tech_sub'),'prespecified 3-test family only')
cat('prev echo HR %:', round(100*mean(A$tachyE),1), '| ICU HR %:', round(100*mean(A$tachyI),1), '\n')
