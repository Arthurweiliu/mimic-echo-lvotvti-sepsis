## 74_lrt_mi.R — non-linearity and overall tests within each imputation
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(splines))
covs <- 'age + sex + sofa + charlson + lactate_filled2 + creat'
z1 <- readRDS('mi_comp_H1.rds'); z2 <- readRDS('mi_comp_H2.rds')

test <- function(z, lin, spl, label, df=2) {
  ps <- sapply(z$comp, function(dk) {
    f0 <- coxph(as.formula(paste('Surv(start,fu_days,event90) ~', covs, '+', lin)), data=dk)
    f1 <- coxph(as.formula(paste('Surv(start,fu_days,event90) ~', covs, '+', spl)), data=dk)
    D <- 2*(f1$loglik[2] - f0$loglik[2])
    pchisq(D, df, lower.tail=FALSE)
  })
  cat(sprintf('%-6s %-28s median p = %.3g (range %.2g-%.2g), %d/%d imputations p<0.05\n',
              label$set, label$what, median(ps), min(ps), max(ps), sum(ps < 0.05), length(ps)))
  invisible(ps)
}
cat('--- non-linearity (spline vs linear on the same scale) ---\n')
test(z1, 'log(vti)',  'v_s1 + v_s2 + v_s3', list(set='H1', what='VTI'))
test(z1, 'lvef',      'l_s1 + l_s2 + l_s3', list(set='H1', what='LVEF'))
test(z2, 'log(ci)',   'e_s1 + e_s2 + e_s3', list(set='H2', what='eCI'))
cat('\n--- overall association (spline vs clinical only) ---\n')
omb <- function(z, spl, label) {
  ps <- sapply(z$comp, function(dk) {
    f0 <- coxph(as.formula(paste('Surv(start,fu_days,event90) ~', covs)), data=dk)
    f1 <- coxph(as.formula(paste('Surv(start,fu_days,event90) ~', covs, '+', spl)), data=dk)
    pchisq(2*(f1$loglik[2] - f0$loglik[2]), 3, lower.tail=FALSE)
  })
  cat(sprintf('%-6s %-28s median p = %.3g (range %.2g-%.2g)\n', label$set, label$what,
              median(ps), min(ps), max(ps)))
}
omb(z1, 'v_s1 + v_s2 + v_s3', list(set='H1', what='VTI'))
omb(z1, 'l_s1 + l_s2 + l_s3', list(set='H1', what='LVEF'))
omb(z2, 'e_s1 + e_s2 + e_s3', list(set='H2', what='eCI'))
