## 69_figures_mi.R — Figures 3, 4, 5 on the multiply imputed cohort
## Fig 3: per-SD estimates (A) and delta AUC (B)
## Fig 4: spline curves for VTI, LVEF and eCI (pooled over imputations, Rubin variance)
## Fig 5: calibration deciles on the imputation-averaged risk
source("config.R")
suppressMessages(library(survival)); suppressMessages(library(splines))

z1 <- readRDS('mi_comp_H1.rds'); z2 <- readRDS('mi_comp_H2.rds')
t3 <- read.csv('Table3_mi.csv'); t4 <- read.csv('Table4_mi_perSD.csv')
covs <- 'age + sex + sofa + charlson + lactate_filled2 + creat'

## ---------------- Fig 3 ----------------
png('Fig3_headtohead_v5.png', width=3600, height=2500, res=600)
par(mfrow=c(1,2), las=1)
par(mar=c(4.8,5.2,3.0,1.2), mgp=c(2.4,0.6,0))
o <- c('VTI','LVEF','eCI')
est <- t4$est[match(o, t4$X)]; lo <- t4$lo[match(o, t4$X)]; hi <- t4$hi[match(o, t4$X)]
plot(NA, xlim=c(min(lo)-0.05, max(hi)+0.05), ylim=c(0.4,3.6), yaxt='n', bty='n',
     xlab='Adjusted hazard ratio per SD', ylab='')
axis(1, cex.axis=0.9); axis(2, at=3:1, labels=o, cex.axis=0.9)
abline(v=1, lty=2, col='grey50')
cols <- c('#2166ac','#4d9221','#b2182b')
for (i in 1:3) {
  segments(lo[i], 4-i, hi[i], 4-i, lwd=2.8, col=cols[i])
  points(est[i], 4-i, pch=19, cex=1.3, col=cols[i])
}
title('A. Per-SD HR (H2, MI)', adj=0, cex.main=0.95)
par(mar=c(5.2,5.2,3.4,2.2))
h1 <- t3[t3$set=='H1', ]; h2 <- t3[t3$set=='H2', ]
labs <- c('VTI','LVEF','VTI+LVEF','VTI','LVEF','eCI','VTI+LVEF','All')
vals <- c(h1$dAUC[h1$model %in% c('VTI','LVEF','VTI_LVEF')],
          h2$dAUC[h2$model %in% c('VTI','LVEF','eCI','VTI_LVEF','all')])
grp <- c(rep('H1',3), rep('H2',5)); cols2 <- ifelse(grp=='H1','#2166ac','#b2182b')
barplot(vals, names.arg=labs, col=cols2, las=2, cex.names=0.72, ylim=c(0, max(vals)*1.15),
        ylab='Delta AUC', border=NA)
abline(h=0.01, lty=2, col='grey30')
title('B. Delta AUC (MI)', adj=0, cex.main=0.9)
dev.off()
cat('Fig3 done\n')

## ---------------- Fig 4: spline curves, pooled over imputations ----------------
curve_pool <- function(z, vname, tag, grid) {
  X <- z$X; comp <- z$comp
  Bfull <- if (vname == 'vti') ns(log(X$vti), df=3) else if (vname == 'lvef') ns(X$lvef, df=3) else ns(log(X$ci), df=3)
  xmed <- if (vname == 'lvef') median(X$lvef) else median(log(X[[vname]]))
  Bmed <- predict(Bfull, xmed)
  Bg <- predict(Bfull, if (vname == 'lvef') grid else log(grid))
  terms <- paste0(tag, '_s', 1:3)
  thetas <- matrix(NA, length(comp), length(grid)); vars <- thetas
  for (k in seq_along(comp)) {
    dk <- comp[[k]]
    f <- coxph(as.formula(paste('Surv(start,fu_days,event90) ~', covs, '+', paste(terms, collapse='+'))), data=dk)
    b <- coef(f)[terms]; V <- vcov(f)[terms, terms]
    for (j in seq_along(grid)) {
      cj <- as.numeric(Bg[j, ] - Bmed)
      thetas[k, j] <- sum(cj * b)
      vars[k, j] <- as.numeric(t(cj) %*% V %*% cj)
    }
  }
  th <- colMeans(thetas); m <- length(comp)
  Tv <- colMeans(vars) + (1 + 1/m) * apply(thetas, 2, stats::var)
  data.frame(x=grid, hr=exp(th), lo=exp(th - 1.96*sqrt(Tv)), hi=exp(th + 1.96*sqrt(Tv)))
}

g_vti <- seq(8, 40, by=0.5)
g_lv <- seq(10, 80, by=1)
g_ci <- seq(1.0, 7.5, by=0.1)
c_vti <- curve_pool(z1, 'vti', 'v', g_vti)
c_lv  <- curve_pool(z1, 'lvef', 'l', g_lv)
c_ci  <- curve_pool(z2, 'ci', 'e', g_ci)
write.csv(c_vti, 'Fig4_curve_vti.csv', row.names=FALSE)
write.csv(c_lv, 'Fig4_curve_lvef.csv', row.names=FALSE)
write.csv(c_ci, 'Fig4_curve_eci.csv', row.names=FALSE)

png('Fig4_spline_curves_v5.png', width=3600, height=1500, res=600)
par(mfrow=c(1,3), las=1, mar=c(4.6,4.8,2.6,0.8), mgp=c(2.4,0.6,0))
draw <- function(d, xlab, main, ref) {
  plot(NA, xlim=range(d$x), ylim=c(0.4, max(2.2, min(3.2, max(d$hi, na.rm=TRUE)))), bty='n',
       xlab=xlab, ylab='Hazard ratio (95% CI)', las=1, log='y')
  polygon(c(d$x, rev(d$x)), c(d$lo, rev(d$hi)), col='grey88', border=NA)
  lines(d$x, d$hr, lwd=2.6, col='#b2182b')
  abline(h=1, lty=2, col='grey40'); abline(v=ref, lty=3, col='grey60')
  title(main, adj=0, cex.main=0.95)
}
draw(c_vti, 'LVOT-VTI, cm', 'A. VTI', median(z1$X$vti))
draw(c_lv, 'LVEF, %', 'B. LVEF', median(z1$X$lvef))
draw(c_ci, expression('Derived CI, L/min/m'^2), 'C. eCI', median(z2$X$ci))
dev.off()
cat('Fig4 done\n')

## ---------------- Fig 5: calibration on averaged risk ----------------
calib_panel <- function(z, model, main, col) {
  p <- z$avg[, model]; y <- z$X$event90
  q <- cut(p, unique(quantile(p, seq(0, 1, 0.1))), include.lowest=TRUE)
  obs <- tapply(y, q, mean); prd <- tapply(p, q, mean); nn <- tapply(y, q, length)
  se <- sqrt(obs*(1-obs)/nn)
  plot(NA, xlim=c(0, max(p)*1.05), ylim=c(0, 1), bty='n', las=1,
       xlab='Predicted 90-day risk', ylab='Observed 90-day risk', main='')
  abline(0, 1, lty=2, col='grey40')
  segments(prd, pmax(obs-1.96*se, 0), prd, pmin(obs+1.96*se, 1), col=col, lwd=1.6)
  points(prd, obs, pch=19, col=col, cex=0.9)
  title(main, adj=0, cex.main=0.95)
}
png('Fig5_calibration_v5.png', width=2800, height=2500, res=600)
par(mfrow=c(2,2), las=1, mar=c(4.4,4.6,2.4,0.9), mgp=c(2.4,0.6,0))
calib_panel(z1, 'clinical', 'A. H1 clinical', '#2166ac')
calib_panel(z1, 'VTI_LVEF', 'B. H1 VTI+LVEF splines', '#b2182b')
calib_panel(z2, 'clinical', 'C. H2 clinical', '#2166ac')
calib_panel(z2, 'all', 'D. H2 all splines', '#b2182b')
dev.off()
cat('Fig5 done\n')
