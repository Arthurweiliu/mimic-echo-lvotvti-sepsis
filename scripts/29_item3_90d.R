# 29_item3_90d.R — Item-3 incremental metrics on the PRIMARY outcome (90-day mortality)
# plus a risk-category NRI (in addition to the category-free NRI).
{ .a <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num <- function(x) suppressWarnings(as.numeric(x))
d <- read.csv('h1_data.csv', stringsAsFactors = FALSE)
for (v in c('age','sofa','charlson','lactate','creat','fu_days','event90','death28','death90','vti')) d[[v]] <- num(d[[v]])
d$sex <- factor(ifelse(d$sex %in% c('M','F'), d$sex, NA))
covs <- c('age','sex','sofa','charlson','lactate','creat'); cov <- paste(covs, collapse = '+')

aucf <- function(p, y) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1*(n1+1)/2)/(n1*n0) }
dca  <- function(p, y, pt) { n <- length(y); (sum(p >= pt & y == 1) - sum(p >= pt & y == 0) * pt/(1-pt))/n }
nri_free <- function(p0, p1, y) { up <- p1>p0; dn <- p1<p0
  (mean(up[y==1]) - mean(dn[y==1])) + (mean(dn[y==0]) - mean(up[y==0])) }
nri_cat <- function(p0, p1, y, br) { c0 <- cut(p0, br, labels = FALSE); c1 <- cut(p1, br, labels = FALSE)
  e <- y == 1; ne <- y == 0
  (mean(c1[e] > c0[e]) - mean(c1[e] < c0[e])) + (mean(c1[ne] < c0[ne]) - mean(c1[ne] > c0[ne])) }
idi <- function(p0, p1, y) (mean(p1[y==1]) - mean(p1[y==0])) - (mean(p0[y==1]) - mean(p0[y==0]))

run <- function(outcome, br, ptmax, tag) {
  B <- d[complete.cases(d[, c(outcome, covs, 'vti')]) & num(d$vti) > 0, ]
  B$y <- as.numeric(B[[outcome]]); B$vti5 <- B$vti/5
  m0 <- glm(as.formula(paste('y~', cov)), data = B, family = binomial)
  m1 <- glm(as.formula(paste('y~', cov, '+ vti5')), data = B, family = binomial)
  p0 <- predict(m0, type = 'response'); p1 <- predict(m1, type = 'response')
  A0 <- aucf(p0, B$y); A1 <- aucf(p1, B$y)
  nf <- nri_free(p0, p1, B$y); nc <- nri_cat(p0, p1, B$y, br); iv <- idi(p0, p1, B$y)
  cat(sprintf('\n=== %s : N=%d events=%d (%.1f%%) ===\n', tag, nrow(B), sum(B$y), 100*mean(B$y)))
  cat(sprintf('  AUC base=%.3f +VTI=%.3f deltaAUC=%+.4f\n', A0, A1, A1-A0))
  cat(sprintf('  NRI(category-free)=%.3f  NRI(risk-category)=%.3f  IDI=%.4f\n', nf, nc, iv))
  tp <- seq(0.02, ptmax, by = 0.01)
  nb <- data.frame(pt = tp, none = 0, all = sapply(tp, function(t) mean(B$y)-(1-mean(B$y))*t/(1-t)),
    base = sapply(tp, function(t) dca(p0, B$y, t)), base_VTI = sapply(tp, function(t) dca(p1, B$y, t)))
  nb$delta <- nb$base_VTI - nb$base
  write.csv(nb, paste0('Item3_DCA_', tag, '.csv'), row.names = FALSE)
  cat('  DCA delta NB: '); for (t in seq(0.1, min(ptmax,0.5), 0.1)) { i <- which.min(abs(nb$pt-t))
    cat(sprintf('pt=%.2f:%+.4f ', t, nb$delta[i])) }; cat('\n')
  png(paste0('Item3_DCA_', tag, '.png'), width = 1200, height = 950, res = 150)
  plot(nb$pt, nb$base_VTI, type='l', lwd=3, col='#1f6fb2', ylim=range(c(0, nb$all, nb$base, nb$base_VTI), na.rm=TRUE),
       xlab='Threshold probability', ylab='Net benefit', main=paste0('DCA — ', tag))
  lines(nb$pt, nb$base, lwd=2, col='#d1495b'); lines(nb$pt, nb$all, lwd=1.2, lty=2, col='grey40')
  abline(h=0, lwd=1.2, col='grey70')
  legend('topright', c('Clinical + LVOT-VTI','Clinical only','Treat all','Treat none'),
         col=c('#1f6fb2','#d1495b','grey40','grey70'), lwd=c(3,2,1.2,1.2), lty=c(1,1,2,1), bty='n'); dev.off()
  # bootstrap CI for the incremental metrics
  set.seed(20261014); Bn <- 500
  M <- matrix(NA, Bn, 5, dimnames=list(NULL, c('dAUC','NRI_free','NRI_cat','IDI','NBd_20')))
  for (i in 1:Bn) { idx <- sample(nrow(B), replace=TRUE); bi <- B[idx,]
    a <- try(glm(as.formula(paste('y~',cov)), data=bi, family=binomial), silent=TRUE)
    b2 <- try(glm(as.formula(paste('y~',cov,'+ vti5')), data=bi, family=binomial), silent=TRUE)
    if (inherits(a,'try-error')||inherits(b2,'try-error')||length(unique(bi$y))<2) next
    q0 <- predict(a,type='response'); q1 <- predict(b2,type='response')
    M[i,'dAUC'] <- aucf(q1,bi$y)-aucf(q0,bi$y); M[i,'NRI_free'] <- nri_free(q0,q1,bi$y)
    M[i,'NRI_cat'] <- nri_cat(q0,q1,bi$y,br); M[i,'IDI'] <- idi(q0,q1,bi$y)
    M[i,'NBd_20'] <- dca(q1,bi$y,0.20)-dca(q0,bi$y,0.20) }
  ci <- t(apply(M, 2, function(x) quantile(x, c(.025,.975), na.rm=TRUE)))
  est <- c(A1-A0, nf, nc, iv, dca(p1,B$y,0.20)-dca(p0,B$y,0.20))
  res <- data.frame(metric=colnames(M), estimate=round(est,4), boot_lo=round(ci[,1],4), boot_hi=round(ci[,2],4))
  write.csv(res, paste0('Item3_incremental_', tag, '.csv'), row.names = FALSE)
  print(res)
}
run('death28', c(-Inf, 0.10, 0.25, Inf), 0.60, '28d')
run('death90', c(-Inf, 0.10, 0.25, Inf), 0.70, '90d')
