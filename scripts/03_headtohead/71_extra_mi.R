## 71_extra_mi.R — supplementary numbers for the MI Results text
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(splines))
covs <- 'age + sex + sofa + charlson + lactate_filled2 + creat'
auc <- function(y, p) { r <- rank(p); n1 <- sum(y==1); n0 <- sum(y==0); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }

avg_risks <- function(z, forms) {
  X <- z$X; comp <- z$comp; n <- nrow(X)
  P <- matrix(0, n, length(forms), dimnames=list(NULL, names(forms)))
  for (nm in names(forms)) for (k in seq_along(comp)) {
    dk <- comp[[k]]
    f <- coxph(as.formula(forms[[nm]]), data=dk)
    bh <- basehaz(f, centered=FALSE)
    H90 <- approx(bh$time, bh$hazard, 90, rule=2)$y
    Hs <- approx(bh$time, bh$hazard, dk$start, rule=2)$y
    lp <- predict(f, newdata=dk, type='lp', reference='zero')
    P[, nm] <- P[, nm] + (1 - exp(-(H90 - Hs)*exp(lp)))/length(comp)
  }
  P
}
z1 <- readRDS('mi_comp_H1.rds'); z2 <- readRDS('mi_comp_H2.rds')

## linear (untransformed) models in H1, for the "linear coefficients understate" sentence
lin <- list(clinical=paste('Surv(start,fu_days,event90) ~', covs),
            VTI_lin=paste('Surv(start,fu_days,event90) ~', covs, '+ vti'),
            LVEF_lin=paste('Surv(start,fu_days,event90) ~', covs, '+ lvef'),
            both_lin=paste('Surv(start,fu_days,event90) ~', covs, '+ vti + lvef'))
P1 <- avg_risks(z1, lin)
A1 <- sapply(colnames(P1), function(v) auc(z1$X$event90, P1[, v]))
cat('H1 linear models (MI):\n'); print(round(A1, 4))
cat('  dAUC linear VTI', round(A1['VTI_lin'] - A1['clinical'], 4),
    '| linear LVEF', round(A1['LVEF_lin'] - A1['clinical'], 4),
    '| both linear', round(A1['both_lin'] - A1['clinical'], 4), '\n')

## correlation between VTI and eCI in H2
r <- with(z2$X, cor(vti, ci, use='complete.obs'))
cat(sprintf('r(VTI, eCI) in H2 = %.3f | n pairs %d\n', r, sum(complete.cases(z2$X[, c('vti','ci')]))))

## calibration detail: mean predicted vs observed, deciles outside the binomial CI
calib_detail <- function(z, model) {
  p <- z$avg[, model]; y <- z$X$event90
  q <- cut(p, unique(quantile(p, seq(0, 1, 0.1))), include.lowest=TRUE)
  obs <- tapply(y, q, mean); prd <- tapply(p, q, mean); nn <- tapply(y, q, length)
  se <- sqrt(obs*(1-obs)/nn)
  out <- sum(abs(obs - prd) > 1.96*se)
  cat(sprintf('  %s: mean predicted %.4f vs observed %.4f | deciles outside CI: %d/10\n',
              model, mean(p), mean(y), out))
}
cat('\nCalibration detail:\n')
for (m in c('clinical','VTI_LVEF')) calib_detail(z1, m)
for (m in c('clinical','all')) calib_detail(z2, m)
