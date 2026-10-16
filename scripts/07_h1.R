# 07_h1.R — H1: 低 LVOT-VTI（每搏量）与脓毒症 90 天死亡
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num <- function(x) suppressWarnings(as.numeric(x))
d <- read.csv('h1_data.csv', stringsAsFactors=FALSE)
for (v in c('age','sofa','sapsii','charlson','lactate','creat','fu_days','event90','death28','death90',
            'vti','vti_diam','lvef','trv','lad','norepi','los_icu','aki_any','aki3','crrt','mbp_min','hr_mean'))
  d[[v]] <- num(d[[v]])
d$sex <- factor(ifelse(d$sex %in% c('M','F'), d$sex, NA)); d$vent <- factor(d$vent)
suppressMessages(library(survival))
auc <- function(p,y){ y<-as.numeric(y); ok<-!is.na(p)&!is.na(y); p<-p[ok]; y<-y[ok]
  n1<-sum(y==1); n0<-sum(y==0); if(!n1||!n0) return(NA_real_)
  r<-rank(p); (sum(r[y==1])-n1*(n1+1)/2)/(n1*n0) }
Cman <- function(lp,tm,ev){ n<-length(lp); num<-0; den<-0
  for(i in 1:(n-1)) for(j in (i+1):n){ ti<-tm[i]; tj<-tm[j]
    if(ev[i]==1 && ti<tj){ den<-den+1; if(lp[i]>lp[j]) num<-num+1 else if(lp[i]==lp[j]) num<-num+0.5 }
    if(ev[j]==1 && tj<ti){ den<-den+1; if(lp[j]>lp[i]) num<-num+1 else if(lp[i]==lp[j]) num<-num+0.5 } }
  if(den==0) return(NA_real_); num/den }
HRr <- function(fit,term){ b<-summary(fit)$coefficients
  sprintf('%.2f (%.2f-%.2f) p=%.2g', exp(b[term,1]), exp(b[term,1]-1.96*b[term,3]), exp(b[term,1]+1.96*b[term,3]), b[term,5]) }
HRr5 <- function(fit,term){ b<-summary(fit)$coefficients   # per +5 cm (coef is per +1 cm)
  sprintf('%.2f (%.2f-%.2f) p=%.2g', exp(b[term,1]*5), exp((b[term,1]-1.96*b[term,3])*5), exp((b[term,1]+1.96*b[term,3])*5), b[term,5]) }
covs <- c('age','sex','sofa','charlson','lactate','creat')
cov <- paste(covs, collapse='+')
v <- d$vti; cat('VTI 可行(=v>0):', sum(!is.na(v)&v>0), '/', nrow(d), sprintf('(%.1f%%)', mean(!is.na(v)&v>0)*100), '\n')
cat('VTI 分布: 中位', median(v,na.rm=TRUE), 'IQR', paste(quantile(v,c(.25,.75),na.rm=TRUE),collapse='-'), 'range', paste(range(v,na.rm=TRUE),collapse='-'), '\n')
cat('分类: <17cm', sum(v<17,na.rm=TRUE), '| 17-20', sum(v>=17&v<=20,na.rm=TRUE), '| >20', sum(v>20,na.rm=TRUE), '\n')
cat('90d死亡率按 VTI: <17=', round(mean(d$event90[v<17],na.rm=TRUE)*100,1), '% 17-20=',
    round(mean(d$event90[v>=17&v<=20],na.rm=TRUE)*100,1), '% >20=', round(mean(d$event90[v>20],na.rm=TRUE)*100,1),'%\n')

base <- complete.cases(d[,c('event90','fu_days',covs,'vti')]) & d$fu_days>0 & d$vti>0
dd <- d[base,]; cat('\n主分析 N=',nrow(dd),' 事件=',sum(dd$event90==1),'\n')
m0 <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov)), data=dd)
cat('  per +5cm VTI:', HRr5(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti')),data=dd),'vti'),
    ' | per +1cm:', HRr(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti')),data=dd),'vti'),'\n')
dd$vcat <- factor(ifelse(dd$vti<17,'lt17', ifelse(dd$vti<=20,'17to20','gt20')), levels=c('gt20','lt17','17to20'))
m2 <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vcat')), data=dd)
cat('  vs >20cm:  <17:', HRr(m2,'vcatlt17'),' | 17-20:', HRr(m2,'vcat17to20'),'\n')
dd$vq <- relevel(cut(dd$vti, quantile(dd$vti, 0:5/5), labels=c('Q1','Q2','Q3','Q4','Q5'), include.lowest=TRUE), ref='Q5')
m3 <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vq')), data=dd)
cat('  五分位(Q5=最高VTI为参照):\n'); for(q in c('Q1','Q2','Q3','Q4')) cat('    ',q,' ',HRr(m3,paste0('vq',q)),'\n')
m4 <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ splines::ns(vti,3)')), data=dd)
m1 <- coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti')), data=dd)
a<-anova(m0,m4,test='LRT'); cat('  样条(3df) LRT p=', format(a[2,ncol(a)],scientific=TRUE,digits=3), '\n')
cat('  C-index: 基线', round(summary(m0)$concordance[1],3), ' +VTI(连续)', round(summary(m1)$concordance[1],3), '\n')

# 28天 logistic + AUC
b2 <- complete.cases(d[,c('death28',covs,'vti')]) & d$vti>0; dd2 <- d[b2,]
l0 <- glm(as.formula(paste('death28~',cov)),data=dd2,family=binomial)
l1 <- glm(as.formula(paste('death28~',cov,'+ vti')),data=dd2,family=binomial)
cat('  28d AUC 基线',round(auc(predict(l0,type='response'),dd2$death28),3),' +VTI',round(auc(predict(l1,type='response'),dd2$death28),3),'\n')

# 次要结局
for (o in c('aki_any','aki3','crrt')) {covs2 <- c('age','sex','sofa','charlson','lactate'); cov2 <- paste(covs2, collapse='+')
cc <- complete.cases(d[,c(o,covs2,'vti')]) & d$vti>0; ddx <- d[cc,]
  f <- glm(as.formula(paste(o,'~',cov2,'+ vti')),data=ddx,family=binomial); b<-summary(f)$coefficients
  cat(sprintf('  %s: VTI每+5cm OR=%.2f (%.2f-%.2f) p=%.2g\n', o, exp(b['vti',1]*5), exp((b['vti',1]-1.96*b['vti',2])*5), exp((b['vti',1]+1.96*b['vti',2])*5), b['vti',4]))
}
# 亚组
for (s in list(c('vent','1'),c('norepi','>0'))) {
  if (s[2]=='1') sub <- dd[dd$vent=='t',] else sub <- dd[dd$norepi>0,]
  if (nrow(sub)>50) { f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti')),data=sub)
    cat(sprintf('  亚组 %s: per+5cm %s (N=%d)\n', s[1], HRr5(f,'vti'), nrow(sub))) }
}
# 时间验证
if (sum(!is.na(dd$yg))>100) {
  yg <- dd$yg
  if (is.factor(yg) || is.character(yg)) { lv <- sort(unique(yg)); cat('  yg 组:', paste(lv,collapse=','),'\n')
    tr <- dd[yg %in% lv[1:max(1,floor(length(lv)/2))],]; va <- dd[!yg %in% lv[1:max(1,floor(length(lv)/2))],]
    if(nrow(va)>50 && sum(va$event90)>20){ f<-coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti')),data=tr)
      lp <- predict(f, newdata=va, type='lp'); cv <- Cman(lp, va$fu_days, va$event90)
      cat('  时间验证 推导C=',round(summary(f)$concordance[1],3),' 验证C=',round(cv,3),' (推导N=',nrow(tr),' 验证N=',nrow(va),')\n') }
  }
} else cat('  (anchor_year_group 缺失，跳过时间验证)\n')
# KM figure is produced by 18_km_tdhr.R (Figure 2, with number-at-risk) — not duplicated here.
