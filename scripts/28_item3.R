# 28_item3.R — Item-3 incremental package (no new data):
#   (1) Landmark (72 h) as co-primary analysis  -> removes immortal-time bias
#   (2) AF-stratified analysis of LVOT-VTI (effect modification, made explicit)
#   (3) Decision-curve analysis (DCA) + NRI + IDI  -> decision-level incremental value
# Cohort and covariates identical to the primary analysis (07_h1.R / 27_validity.R).
{ .a <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }

num <- function(x) suppressWarnings(as.numeric(x))
d <- read.csv('h1_data.csv', stringsAsFactors = FALSE)
for (v in c('age','sofa','charlson','lactate','creat','fu_days','event90','death28','death90','vti'))
  d[[v]] <- num(d[[v]])
d$sex <- factor(ifelse(d$sex %in% c('M','F'), d$sex, NA))
suppressMessages(library(survival))

## ---- AF flag (identical to 27_validity.R) ----
af1 <- read.csv('af_icu.csv', stringsAsFactors = FALSE)
af2 <- read.csv('af_icd.csv', stringsAsFactors = FALSE)
d$stay_id <- as.character(d$stay_id); af1$stay_id <- as.character(af1$stay_id); af2$stay_id <- as.character(af2$stay_id)
d <- merge(d, af1[, c('stay_id','af_icu')], by = 'stay_id', all.x = TRUE)
d <- merge(d, af2[, c('stay_id','af_icd')], by = 'stay_id', all.x = TRUE)
d$af <- as.numeric((!is.na(d$af_icu) & d$af_icu == 1) | (!is.na(d$af_icd) & d$af_icd == 1))

covs <- c('age','sex','sofa','charlson','lactate','creat'); cov <- paste(covs, collapse = '+')
ccs <- c('event90','fu_days', covs, 'vti')
A <- d[complete.cases(d[, ccs]) & num(d$fu_days) > 0 & num(d$vti) > 0, ]
A$vlow <- as.numeric(A$vti < 17)
cat('=== Primary analysis set ===\n'); cat('N =', nrow(A), ' events(90d) =', sum(A$event90 == 1),
    ' AF =', sum(A$af), sprintf('(%.1f%%)\n', 100 * mean(A$af)), '\n\n')

hrline <- function(x, term, lab, per5 = FALSE) {
  f <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+', term)), data = x)
  b <- summary(f)$coefficients
  k <- if (per5) 5 else 1
  cat(sprintf('%-52s N=%5d  HR=%.2f (%.2f-%.2f) p=%.2g\n', lab, nrow(x),
    exp(b[nrow(b), 1] * k), exp((b[nrow(b), 1] - 1.96 * b[nrow(b), 3]) * k),
    exp((b[nrow(b), 1] + 1.96 * b[nrow(b), 3]) * k), b[nrow(b), 5]))
  invisible(summary(f))
}

## ======================= (1) LANDMARK AS CO-PRIMARY =========================
cat('=== (1) Landmark analysis (time origin = L days after ICU admission) ===\n')
lm_tab <- list()
for (L in c(0, 1, 2, 3, 7)) {
  Al <- A[A$fu_days > L, ]
  Al$tl <- Al$fu_days - L
  f <- coxph(as.formula(paste('Surv(tl,event90)~', cov, '+ vlow')), data = Al)
  b <- summary(f)$coefficients
  f2 <- coxph(as.formula(paste('Surv(tl,event90)~', cov, '+ vti')), data = Al)
  b2 <- summary(f2)$coefficients
  lm_tab[[length(lm_tab) + 1]] <- data.frame(Landmark_days = L, N = nrow(Al), Events = sum(Al$event90 == 1),
    HR_vlow = round(exp(b['vlow', 1]), 2), lo = round(exp(b['vlow', 1] - 1.96 * b['vlow', 3]), 2),
    hi = round(exp(b['vlow', 1] + 1.96 * b['vlow', 3]), 2), p = format.pval(b['vlow', 5], digits = 2),
    HR_per5 = round(exp(b2['vti', 1] * 5), 3), lo5 = round(exp((b2['vti', 1] - 1.96 * b2['vti', 3]) * 5), 3),
    hi5 = round(exp((b2['vti', 1] + 1.96 * b2['vti', 3]) * 5), 3))
  cat(sprintf('  L=%dd  N=%5d ev=%4d  vlow HR=%.2f (%.2f-%.2f) p=%.2g | per+5cm HR=%.3f (%.3f-%.3f)\n',
    L, nrow(Al), sum(Al$event90 == 1), exp(b['vlow',1]), exp(b['vlow',1]-1.96*b['vlow',3]),
    exp(b['vlow',1]+1.96*b['vlow',3]), b['vlow',5], exp(b2['vti',1]*5),
    exp((b2['vti',1]-1.96*b2['vti',3])*5), exp((b2['vti',1]+1.96*b2['vti',3])*5)))
}
lm_tab <- do.call(rbind, lm_tab); write.csv(lm_tab, 'Item3_landmark.csv', row.names = FALSE)

## ======================= (2) AF STRATIFICATION ==============================
cat('\n=== (2) Effect modification by atrial fibrillation (AF) ===\n')
f <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ vlow*af')), data = A)
b <- summary(f)$coefficients
it <- grep('^vlow:af$', rownames(b))
cat(sprintf('  Interaction vlow x AF: HR=%.2f (%.2f-%.2f) p=%.3f\n',
  exp(b[it,1]), exp(b[it,1]-1.96*b[it,3]), exp(b[it,1]+1.96*b[it,3]), b[it,5]))
af_tab <- list()
for (g in c(0, 1)) {
  s <- A[A$af == g, ]
  ff <- coxph(as.formula(paste('Surv(fu_days,event90)~', cov, '+ vlow')), data = s)
  bb <- summary(ff)$coefficients
  lab <- if (g == 0) 'No AF' else 'AF'
  cat(sprintf('  %-6s N=%5d ev=%4d  vlow HR=%.2f (%.2f-%.2f)\n', lab, nrow(s), sum(s$event90==1),
    exp(bb['vlow',1]), exp(bb['vlow',1]-1.96*bb['vlow',3]), exp(bb['vlow',1]+1.96*bb['vlow',3])))
  af_tab[[length(af_tab)+1]] <- data.frame(Stratum = lab, N = nrow(s), Events = sum(s$event90==1),
    HR = round(exp(bb['vlow',1]),2), lo = round(exp(bb['vlow',1]-1.96*bb['vlow',3]),2),
    hi = round(exp(bb['vlow',1]+1.96*bb['vlow',3]),2))
}
write.csv(do.call(rbind, af_tab), 'Item3_af.csv', row.names = FALSE)

## ======================= (3) DCA + NRI/IDI =================================
cat('\n=== (3) Decision-level incremental value (28-day mortality) ===\n')
# binary 28-day mortality; complete cases as in the primary model
B <- d[complete.cases(d[, c('death28', covs, 'vti')]) & num(d$vti) > 0, ]
B$y <- as.numeric(B$death28)
B$vti5 <- B$vti / 5
m0 <- glm(as.formula(paste('y~', cov)), data = B, family = binomial)
m1 <- glm(as.formula(paste('y~', cov, '+ vti5')), data = B, family = binomial)
p0 <- predict(m0, type = 'response'); p1 <- predict(m1, type = 'response')
cat('  N =', nrow(B), ' events =', sum(B$y), sprintf('(%.1f%%)\n', 100*mean(B$y)))

aucf <- function(p, y) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
A0 <- aucf(p0, B$y); A1 <- aucf(p1, B$y)
cat(sprintf('  AUC base=%.3f  +VTI=%.3f  deltaAUC=%.4f\n', A0, A1, A1 - A0))

# --- DCA net benefit ---
dca <- function(p, y, pt) { n <- length(y)
  (sum(p >= pt & y == 1) - sum(p >= pt & y == 0) * pt / (1 - pt)) / n }
tp <- seq(0.02, 0.60, by = 0.01)
nb <- data.frame(pt = tp,
  none = 0,
  all = sapply(tp, function(t) mean(B$y) - (1 - mean(B$y)) * t / (1 - t)),
  base = sapply(tp, function(t) dca(p0, B$y, t)),
  base_VTI = sapply(tp, function(t) dca(p1, B$y, t)))
nb$delta <- nb$base_VTI - nb$base
write.csv(nb, 'Item3_DCA.csv', row.names = FALSE)
cat('  DCA net benefit (28d) at selected thresholds:\n')
for (t in c(0.10, 0.15, 0.20, 0.25, 0.30, 0.40)) {
  i <- which.min(abs(nb$pt - t))
  cat(sprintf('    pt=%.2f  base=%.4f  +VTI=%.4f  delta=%+.4f  (all=%.4f)\n', t, nb$base[i], nb$base_VTI[i], nb$delta[i], nb$all[i]))
}
png('Item3_DCA_28d.png', width = 1200, height = 950, res = 150)
plot(nb$pt, nb$base_VTI, type = 'l', lwd = 3, col = '#1f6fb2', ylim = range(c(0, nb$all, nb$base, nb$base_VTI), na.rm = TRUE),
     xlab = 'Threshold probability', ylab = 'Net benefit', main = 'Decision curve analysis (28-day mortality)')
lines(nb$pt, nb$base, lwd = 2, col = '#d1495b')
lines(nb$pt, nb$all, lwd = 1.2, lty = 2, col = 'grey40')
abline(h = 0, lwd = 1.2, col = 'grey70')
legend('topright', c('Clinical + LVOT-VTI','Clinical only','Treat all','Treat none'),
       col = c('#1f6fb2','#d1495b','grey40','grey70'), lwd = c(3,2,1.2,1.2), lty = c(1,1,2,1), bty = 'n')
dev.off()

# --- NRI (category-free) + IDI ---
nri_idi <- function(p0, p1, y) {
  up <- p1 > p0; dn <- p1 < p0
  nri <- (mean(up[y == 1]) - mean(dn[y == 1])) + (mean(dn[y == 0]) - mean(up[y == 0]))
  idi <- (mean(p1[y == 1]) - mean(p1[y == 0])) - (mean(p0[y == 1]) - mean(p0[y == 0]))
  c(NRI = nri, IDI = idi)
}
pt_est <- nri_idi(p0, p1, B$y)
cat(sprintf('  Category-free NRI = %.3f ; IDI = %.4f\n', pt_est['NRI'], pt_est['IDI']))

# bootstrap CI
set.seed(20261014); Bn <- 500
boot <- matrix(NA, Bn, 5, dimnames = list(NULL, c('dAUC','NRI','IDI','NBdelta_20','NBdelta_30')))
for (i in 1:Bn) {
  idx <- sample(nrow(B), replace = TRUE); bi <- B[idx, ]
  a <- try(glm(as.formula(paste('y~', cov)), data = bi, family = binomial), silent = TRUE)
  b2 <- try(glm(as.formula(paste('y~', cov, '+ vti5')), data = bi, family = binomial), silent = TRUE)
  if (inherits(a,'try-error') || inherits(b2,'try-error') || length(unique(bi$y)) < 2) next
  q0 <- predict(a, type = 'response'); q1 <- predict(b2, type = 'response')
  boot[i,'dAUC'] <- aucf(q1, bi$y) - aucf(q0, bi$y)
  ni <- nri_idi(q0, q1, bi$y); boot[i,'NRI'] <- ni['NRI']; boot[i,'IDI'] <- ni['IDI']
  boot[i,'NBdelta_20'] <- dca(q1, bi$y, 0.20) - dca(q0, bi$y, 0.20)
  boot[i,'NBdelta_30'] <- dca(q1, bi$y, 0.30) - dca(q0, bi$y, 0.30)
}
ci <- t(apply(boot, 2, function(x) quantile(x, c(.025, .975), na.rm = TRUE)))
inc <- data.frame(metric = c('deltaAUC','NRI (category-free)','IDI','delta NB @pt=0.20','delta NB @pt=0.30'),
  estimate = round(c(A1 - A0, pt_est['NRI'], pt_est['IDI'],
    dca(p1,B$y,0.20)-dca(p0,B$y,0.20), dca(p1,B$y,0.30)-dca(p0,B$y,0.30)), 4),
  boot_lo = round(ci[,1], 4), boot_hi = round(ci[,2], 4))
write.csv(inc, 'Item3_incremental.csv', row.names = FALSE)
cat('  Bootstrap 95% CI (500x):\n'); print(inc)

sink('Item3_results.txt')
cat('Item-3 incremental package results\nGenerated:', format(Sys.time()), '\n\n')
cat('Primary set N =', nrow(A), ' events =', sum(A$event90 == 1), ' AF =', sum(A$af), '\n\n')
cat('--- Landmark ---\n'); print(lm_tab)
cat('\n--- AF strata ---\n'); print(do.call(rbind, af_tab))
cat('\n--- Incremental value (28d) ---\n'); print(inc)
sink()
cat('\nSaved: Item3_landmark.csv, Item3_af.csv, Item3_DCA.csv, Item3_DCA_28d.png, Item3_incremental.csv, Item3_results.txt\n')
