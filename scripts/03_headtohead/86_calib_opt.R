## 86_calib_opt.R — optimism-corrected calibration slope (one set per run)
## Usage: $env:SET='H1'; Rscript 86_calib_opt.R
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(splines))
SET <- Sys.getenv('SET', 'H1'); B <- as.integer(Sys.getenv('BOOT_REPS', '200'))
set.seed(20261105)
covs <- 'age + sex + sofa + charlson + lactate_filled2 + creat'
form <- if (SET == 'H1')
  paste('Surv(start,fu_days,event90) ~', covs, '+ v_s1 + v_s2 + v_s3 + l_s1 + l_s2 + l_s3') else
  paste('Surv(start,fu_days,event90) ~', covs,
        '+ v_s1 + v_s2 + v_s3 + l_s1 + l_s2 + l_s3 + e_s1 + e_s2 + e_s3')

z <- readRDS(if (SET == 'H1') 'mi_comp_H1.rds' else 'mi_comp_H2.rds')
X <- z$X; comp <- z$comp; n <- nrow(X); y <- X$event90

pred <- function(f, dk) {
  bh <- basehaz(f, centered=FALSE)
  H90 <- approx(bh$time, bh$hazard, 90, rule=2)$y
  Hs <- approx(bh$time, bh$hazard, dk$start, rule=2)$y
  lp <- predict(f, newdata=dk, type='lp', reference='zero')
  pmin(pmax(1 - exp(-(H90 - Hs)*exp(lp)), 1e-6), 1 - 1e-6)
}
slope <- function(yy, p) unname(coef(glm(yy ~ qlogis(p), family=binomial))[2])

## apparent slope on the imputation-averaged risk
P <- rowMeans(sapply(comp, function(dk) { f <- coxph(as.formula(form), data=dk); pred(f, dk) }))
app <- slope(y, P)
cat(sprintf('%s apparent calibration slope = %.3f\n', SET, app))

opt <- numeric(B); t0 <- Sys.time()
for (b in seq_len(B)) {
  idx <- sample.int(n, n, replace=TRUE)
  while (sum(y[idx]) < 20) idx <- sample.int(n, n, replace=TRUE)
  sb <- sf <- numeric(length(comp))
  for (k in seq_along(comp)) {
    dk <- comp[[k]][idx, , drop=FALSE]
    f <- coxph(as.formula(form), data=dk)
    sb[k] <- slope(dk$event90, pred(f, dk))
    sf[k] <- slope(comp[[k]]$event90, pred(f, comp[[k]]))
  }
  opt[b] <- mean(sb) - mean(sf)
  if (b %% 20 == 0) cat(sprintf('  %d/%d (%.1f min)\n', b, B,
                                as.numeric(difftime(Sys.time(), t0, units='mins'))))
}
cat(sprintf('%s mean optimism = %.3f | corrected slope = %.3f\n', SET, mean(opt), app - mean(opt)))
write.csv(data.frame(set=SET, apparent=app, optimism=mean(opt), corrected=app - mean(opt)),
          paste0('TableS_calib_opt_', SET, '.csv'), row.names=FALSE)
