## 73_optimism.R — bootstrap optimism correction for the combined models (one set per run)
## Usage: $env:SET='H1'; Rscript 73_optimism.R
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
X <- z$X; comp <- z$comp; n <- nrow(X)
cc <- function(d, lp) as.numeric(concordance(Surv(d$start, d$fu_days, d$event90) ~ lp, reverse=TRUE)$concordance)

## apparent C on the full data (average over imputations)
app <- mean(sapply(comp, function(dk) cc(dk, predict(coxph(as.formula(form), data=dk), newdata=dk, type='lp'))))
cat(sprintf('%s apparent C = %.4f\n', SET, app))

opt <- numeric(B); t0 <- Sys.time()
for (b in seq_len(B)) {
  idx <- sample.int(n, n, replace=TRUE)
  while (sum(X$event90[idx]) < 20) idx <- sample.int(n, n, replace=TRUE)
  cb <- cf <- numeric(length(comp))
  for (k in seq_along(comp)) {
    dk <- comp[[k]][idx, , drop=FALSE]
    f <- coxph(as.formula(form), data=dk)
    cb[k] <- cc(dk, predict(f, newdata=dk, type='lp'))
    cf[k] <- cc(comp[[k]], predict(f, newdata=comp[[k]], type='lp'))
  }
  opt[b] <- mean(cb) - mean(cf)
  if (b %% 20 == 0) cat(sprintf('  %d/%d (%.1f min)\n', b, B,
                                as.numeric(difftime(Sys.time(), t0, units='mins'))))
}
cat(sprintf('%s mean optimism = %.4f\n', SET, mean(opt)))

cat(sprintf('%s optimism-corrected C = %.4f\n', SET, app - mean(opt)))
write.csv(data.frame(set=SET, apparent=app, optimism=mean(opt), corrected=app - mean(opt)),
          paste0('TableS_optimism_', SET, '.csv'), row.names=FALSE)
