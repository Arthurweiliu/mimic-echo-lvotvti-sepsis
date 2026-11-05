## 68_mi_delta.R — incremental value on the multiply imputed cohort (one set per run)
## Usage: $env:SET='H1'; $env:BOOT_REPS='200'; Rscript 68_mi_delta.R
## Paired bootstrap: patients resampled, both models refitted in all 20 completed datasets.
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(splines))
SET <- Sys.getenv('SET', 'H1'); B <- as.integer(Sys.getenv('BOOT_REPS', '200'))
set.seed(20261105)

covs <- 'age + sex + sofa + charlson + lactate_filled2 + creat'
MOD <- list(
  clinical = paste('Surv(start,fu_days,event90) ~', covs),
  new = if (SET == 'H1') paste('Surv(start,fu_days,event90) ~', covs, '+ v_s1 + v_s2 + v_s3 + l_s1 + l_s2 + l_s3')
        else paste('Surv(start,fu_days,event90) ~', covs,
                   '+ v_s1 + v_s2 + v_s3 + l_s1 + l_s2 + l_s3 + e_s1 + e_s2 + e_s3'))

auc <- function(y, p) { r <- rank(p); n1 <- sum(y==1); n0 <- sum(y==0); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }
nri_free <- function(y, po, pn) {
  ev <- y == 1; ne <- !ev
  (mean(pn[ev] > po[ev]) - mean(pn[ev] < po[ev])) + (mean(pn[ne] < po[ne]) - mean(pn[ne] > po[ne]))
}
nri_cat <- function(y, po, pn, cuts=c(0.10, 0.25)) {
  co <- cut(po, c(-Inf, cuts, Inf)); cn <- cut(pn, c(-Inf, cuts, Inf))
  up <- as.integer(cn) > as.integer(co); dn <- as.integer(cn) < as.integer(co)
  ev <- y == 1; ne <- !ev
  (mean(up[ev]) - mean(dn[ev])) + (mean(dn[ne]) - mean(up[ne]))
}
idi <- function(y, po, pn) {
  ev <- y == 1; ne <- !ev
  (mean(pn[ev]) - mean(pn[ne])) - (mean(po[ev]) - mean(po[ne]))
}
nb <- function(y, p, pt) { tp <- mean(p >= pt & y == 1); fp <- mean(p >= pt & y == 0); tp - fp*(pt/(1-pt)) }

z <- readRDS(if (SET == 'H1') 'mi_comp_H1.rds' else 'mi_comp_H2.rds')
X <- z$X; comp <- z$comp; n <- nrow(X); y <- X$event90
cat(sprintf('%s: n = %d, deaths = %d, imputations = %d, bootstrap reps = %d\n', SET, n, sum(y), length(comp), B))

point <- function(idx) {
  P <- matrix(0, length(idx), 2, dimnames=list(NULL, c('clinical','new')))
  for (nm in c('clinical','new')) {
    for (k in seq_along(comp)) {
      dk <- comp[[k]][idx, , drop=FALSE]
      f <- coxph(as.formula(MOD[[nm]]), data=dk)
      bh <- basehaz(f, centered=FALSE)
      H90 <- approx(bh$time, bh$hazard, 90, rule=2)$y
      Hs <- approx(bh$time, bh$hazard, dk$start, rule=2)$y
      lp <- predict(f, newdata=dk, type='lp', reference='zero')
      P[, nm] <- P[, nm] + (1 - exp(-(H90 - Hs) * exp(lp))) / length(comp)
    }
  }
  P
}

P0 <- point(seq_len(n))
A0 <- c(clinical=auc(y, P0[,'clinical']), new=auc(y, P0[,'new']))
d0 <- c(dAUC      = unname(A0['new'] - A0['clinical']),
        NRI_free  = unname(nri_free(y, P0[,'clinical'], P0[,'new'])),
        NRI_cat   = unname(nri_cat(y, P0[,'clinical'], P0[,'new'])),
        IDI       = unname(idi(y, P0[,'clinical'], P0[,'new'])))
cat(sprintf('AUC %.4f -> %.4f | dAUC %.4f | NRI-free %.3f | NRI-cat %.3f | IDI %.4f\n',
            A0['clinical'], A0['new'], d0['dAUC'], d0['NRI_free'], d0['NRI_cat'], d0['IDI']))
ths <- seq(0.20, 0.50, by=0.05)
dca0 <- data.frame(threshold=ths, clinical=sapply(ths, function(t) nb(y, P0[,'clinical'], t)),
                   new=sapply(ths, function(t) nb(y, P0[,'new'], t)))
print(round(dca0, 4))

res <- matrix(NA_real_, B, 4, dimnames=list(NULL, names(d0)))
t0 <- Sys.time()
for (b in seq_len(B)) {
  idx <- sample.int(n, n, replace=TRUE)
  while (sum(y[idx]) < 20) idx <- sample.int(n, n, replace=TRUE)
  Pb <- point(idx); yb <- y[idx]
  res[b, ] <- c(auc(yb, Pb[,'new']) - auc(yb, Pb[,'clinical']),
                nri_free(yb, Pb[,'clinical'], Pb[,'new']),
                nri_cat(yb, Pb[,'clinical'], Pb[,'new']),
                idi(yb, Pb[,'clinical'], Pb[,'new']))
  if (b %% 10 == 0) cat(sprintf('  %d/%d (%.1f min)\n', b, B,
                                as.numeric(difftime(Sys.time(), t0, units='mins'))))
}
ci <- apply(res, 2, quantile, c(.025, .975), na.rm=TRUE)
out <- data.frame(set=SET, metric=names(d0), est=unname(d0), lo=unname(ci[1, ]), hi=unname(ci[2, ]))
print(round(out[, -c(1,2)], 4))
write.csv(out, paste0('TableS_incremental_mi_', SET, '.csv'), row.names=FALSE)
write.csv(dca0, paste0('TableS_dca_mi_', SET, '.csv'), row.names=FALSE)
cat('\nwritten', paste0('TableS_incremental_mi_', SET, '.csv'), '\n')
