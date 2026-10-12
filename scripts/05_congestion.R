# 05_congestion.R — 右心/静脉淤血负担评分 vs AKI/RRT/死亡（base R + survival）
setwd('F:/openclaw_workspace/sepsis_echo_lit')
d <- read.csv('congestion_analysis.csv', stringsAsFactors = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))
for (v in c('age','sofa','sapsii','charlson','lactate','creat','days','death28','death90',
            'aki_any','aki3','crrt','score3','score4','tr_velocity','la_dimen','la_vol','rvi_bsa','rv_diam','lvef','norepi'))
  d[[v]] <- num(d[[v]])
d$sex <- factor(ifelse(d$sex %in% c('M','F'), d$sex, NA))
cat('N total =', nrow(d), '\n')
cat('可用性: TR', sum(!is.na(d$tr_velocity)), '| LA径', sum(!is.na(d$la_dimen)),
    '| LAVI', sum(!is.na(d$rvi_bsa)), '| RV径', sum(!is.na(d$rv_diam)), '| score3 非缺失', sum(!is.na(d$score3)), '\n\n')

auc <- function(p, y){ y <- as.numeric(y); ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]
  n1 <- sum(y==1); n0 <- sum(y==0); if(n1==0||n0==0) return(NA_real_)
  r <- rank(p); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }

cat('== 各评分档的事件率(%) ==\n')
for (o in c('aki_any','aki3','crrt','death28','death90')) {
  cat(sprintf('%-9s', o)); r <- tapply(d[[o]], d$score3, function(x) round(mean(x, na.rm=TRUE)*100,1))
  cat(paste(names(r), r, sep='=', collapse='  '), '\n')
}

base_cov <- c('age','sex','sofa','charlson','lactate','creat')
run_logit <- function(outcome) {
  cc <- complete.cases(d[, c(outcome, base_cov, 'score3')])
  dd <- d[cc, ]
  m0 <- glm(as.formula(paste(outcome, '~', paste(base_cov, collapse='+'))), data=dd, family=binomial)
  m1 <- glm(as.formula(paste(outcome, '~', paste(base_cov, collapse='+'), '+ score3')), data=dd, family=binomial)
  cat(sprintf('\n== %s (完整病例 N=%d, 事件=%d) ==\n', outcome, nrow(dd), sum(dd[[outcome]]==1)))
  b <- summary(m1)$coefficients
  cat(sprintf('score3 每+1: OR=%.2f (%.2f-%.2f), p=%.3g\n', exp(b['score3',1]), exp(b['score3',1]-1.96*b['score3',2]), exp(b['score3',1]+1.96*b['score3',2]), b['score3',4]))
  cat(sprintf('AUC 基线=%.3f  基线+score3=%.3f  增量=%.3f\n', auc(predict(m0,type='response'), dd[[outcome]]), auc(predict(m1,type='response'), dd[[outcome]]), auc(predict(m1,type='response'), dd[[outcome]])-auc(predict(m0,type='response'), dd[[outcome]])))
  # 分类
  dd$s3f <- factor(dd$score3)
  m2 <- glm(as.formula(paste(outcome, '~', paste(base_cov, collapse='+'), '+ s3f')), data=dd, family=binomial)
  bb <- summary(m2)$coefficients; rn <- grep('^s3f', rownames(bb), value=TRUE)
  for (r in rn) cat(sprintf('  %s vs 0: OR=%.2f (%.2f-%.2f) p=%.3g\n', r, exp(bb[r,1]), exp(bb[r,1]-1.96*bb[r,2]), exp(bb[r,1]+1.96*bb[r,2]), bb[r,4]))
}
for (o in c('aki_any','aki3','crrt')) run_logit(o)

# Cox: 90天死亡
suppressMessages(library(survival))
d$fu_days <- num(d$fu_days); d$event90 <- num(d$event90)
cc <- complete.cases(d[, c('event90','fu_days',base_cov,'score3')]) & d$fu_days > 0
dd <- d[cc, ]
c0 <- coxph(Surv(fu_days, event90) ~ age+sex+sofa+charlson+lactate+creat, data=dd)
c1 <- coxph(Surv(fu_days, event90) ~ age+sex+sofa+charlson+lactate+creat+score3, data=dd)
cat(sprintf('\n== Cox 90天死亡 (N=%d, 事件=%d, 90d死亡%%=%.1f) ==\n', nrow(dd), sum(dd$event90==1), mean(dd$event90)*100))
b <- summary(c1)$coefficients
cat(sprintf('score3 每+1: HR=%.2f (%.2f-%.2f), p=%.3g\n', exp(b['score3',1]), exp(b['score3',1]-1.96*b['score3',3]), exp(b['score3',1]+1.96*b['score3',3]), b['score3',5]))
cat(sprintf('C-index 基线=%.3f  基线+score3=%.3f\n', summary(c0)$concordance[1], summary(c1)$concordance[1]))
# 三组分单独
c2 <- coxph(Surv(fu_days, event90) ~ age+sex+sofa+charlson+lactate+creat+highTR+LAabn+RVabn, data=dd)
b2 <- summary(c2)$coefficients
cat('\n-- 组分分别(校正) --\n')
for (r in c('highTR','LAabn','RVabn')) if (r %in% rownames(b2)) cat(sprintf('%-7s HR=%.2f (%.2f-%.2f) p=%.3g\n', r, exp(b2[r,1]), exp(b2[r,1]-1.96*b2[r,3]), exp(b2[r,1]+1.96*b2[r,3]), b2[r,5]))
