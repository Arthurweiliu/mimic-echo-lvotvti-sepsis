## 75_incremental_lrt.R — does adding each marker to the other improve fit (per imputation)
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(splines))
covs <- 'age + sex + sofa + charlson + lactate_filled2 + creat'
z1 <- readRDS('mi_comp_H1.rds'); z2 <- readRDS('mi_comp_H2.rds')
cmp <- function(z, base_terms, add_terms, label, df) {
  ps <- sapply(z$comp, function(dk) {
    f0 <- coxph(as.formula(paste('Surv(start,fu_days,event90) ~', covs, '+', base_terms)), data=dk)
    f1 <- coxph(as.formula(paste('Surv(start,fu_days,event90) ~', covs, '+', base_terms, '+', add_terms)), data=dk)
    pchisq(2*(f1$loglik[2] - f0$loglik[2]), df, lower.tail=FALSE)
  })
  cat(sprintf('%-6s %-46s median p = %.3g (max %.2g)\n', label$set, label$what, median(ps), max(ps)))
}
vs <- 'v_s1 + v_s2 + v_s3'; ls <- 'l_s1 + l_s2 + l_s3'; es <- 'e_s1 + e_s2 + e_s3'
cat('--- incremental value within H1 (3 df each) ---\n')
cmp(z1, vs, ls, list(set='H1', what='LVEF added to VTI'), 3)
cmp(z1, ls, vs, list(set='H1', what='VTI added to LVEF'), 3)
cat('--- incremental value within H2 (3 df each) ---\n')
cmp(z2, paste(vs, '+', ls), es, list(set='H2', what='eCI added to VTI + LVEF'), 3)
cmp(z2, paste(vs, '+', es), ls, list(set='H2', what='LVEF added to VTI + eCI'), 3)
cmp(z2, paste(ls, '+', es), vs, list(set='H2', what='VTI added to LVEF + eCI'), 3)
