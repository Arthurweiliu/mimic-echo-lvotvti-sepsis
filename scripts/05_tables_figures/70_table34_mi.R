## 70_table34_mi.R — assemble Table 3 and Table 4 for the MI cohort
source("config.R")
t3 <- read.csv('Table3_mi.csv', stringsAsFactors=FALSE)
t4 <- read.csv('Table4_mi_perSD.csv', stringsAsFactors=FALSE, row.names=1)

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

add_inc <- function(z, tab, set) {
  X <- z$X; y <- X$event90; avg <- z$avg
  ref <- avg[, 'clinical']
  tab$NRI_free <- NA; tab$NRI_cat <- NA; tab$IDI <- NA
  for (i in seq_len(nrow(tab))) {
    nm <- tab$model[i]; if (nm == 'clinical') next
    tab$NRI_free[i] <- nri_free(y, ref, avg[, nm])
    tab$NRI_cat[i]  <- nri_cat(y, ref, avg[, nm])
    tab$IDI[i]      <- idi(y, ref, avg[, nm])
  }
  tab
}
t3H1 <- add_inc(readRDS('mi_comp_H1.rds'), t3[t3$set=='H1', ], 'H1')
t3H2 <- add_inc(readRDS('mi_comp_H2.rds'), t3[t3$set=='H2', ], 'H2')
tab <- rbind(t3H1, t3H2)
write.csv(tab, 'Table3_mi_final.csv', row.names=FALSE)

lab <- c(clinical='Clinical', VTI='VTI spline', LVEF='LVEF spline', eCI='eCI spline',
         VTI_LVEF='VTI + LVEF splines', VTI_eCI='VTI + eCI splines',
         LVEF_eCI='LVEF + eCI splines', all='All splines')
fmt <- function(x, d=4) formatC(x, format='f', digits=d)
md <- c('### Table 3. Head-to-head models in the multiply imputed cohort (natural cubic splines, 3 df)', '',
        '| Model | C-index | AUC | delta AUC | Calibration slope | Intercept | Brier | NRI-cat | IDI |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|')
for (s in c('H1','H2')) {
  tt <- tab[tab$set==s, ]
  md <- c(md, sprintf('| **%s (n = %s; %s deaths)** | | | | | | | | |', s,
                      format(tt$n[1], big.mark=','), format(tt$deaths[1], big.mark=',')))
  for (i in seq_len(nrow(tt))) {
    r <- tt[i, ]
    md <- c(md, sprintf('| %s | %s | %s | %s | %s | %s | %s | %s | %s |',
                        lab[r$model], fmt(r$Cindex), fmt(r$AUC),
                        if (r$model=='clinical') '-' else sprintf('%+.4f', r$dAUC),
                        fmt(r$slope, 3), fmt(r$intercept, 3), fmt(r$Brier),
                        if (is.na(r$NRI_cat)) '-' else fmt(r$NRI_cat, 3),
                        if (is.na(r$IDI)) '-' else fmt(r$IDI, 4)))
  }
}
md <- c(md, '',
  '*Predicted 90-day risks are the mean of the risks from the 20 imputations; discrimination and calibration are computed on that averaged risk. AUC, area under the receiver-operating-characteristic curve for 90-day mortality; C-index, Harrell concordance with left truncation; Brier score at 90 days; calibration slope and intercept from a logistic model of the outcome on the logit of the predicted risk; NRI-cat, three-level net reclassification improvement (<10%, 10-25%, >=25%) versus the clinical model; IDI, integrated discrimination improvement versus the clinical model. The clinical model contains age, sex, SOFA, Charlson index, lactate and creatinine.*')
writeLines(md, 'Table3_mi.md')
cat(md, sep='\n')

md4 <- c('### Table 4. Per-SD effects in the mutually adjusted H2 linear model (MI)', '',
         '| Marker | SD | HR per SD (95% CI) | p |', '|---|---:|---|---:|')
sdl <- c(VTI=log(NA), LVEF=NA, eCI=NA)
lab4 <- c(VTI='VTI', LVEF='LVEF', eCI='eCI')
unit <- c(VTI='0.332 (log cm)', LVEF='15.99%', eCI='0.367 (log L/min/m2)')
for (v in c('VTI','LVEF','eCI')) {
  r <- t4[v, ]
  md4 <- c(md4, sprintf('| %s | %s | %s (%s-%s) | %s |', lab4[v], unit[v],
                        fmt(r$est, 3), fmt(r$lo, 3), fmt(r$hi, 3),
                        formatC(r$p, format='g', digits=2)))
}
md4 <- c(md4, '', '*Pooled across 20 imputations (Rubin variance). SD is the standard deviation of the modelled term (log-transformed for VTI and eCI, linear for LVEF): 0.332 log units, 15.99%, and 0.367 log units, respectively.*')
writeLines(md4, 'Table4_mi.md')
cat('\n', md4, sep='\n')
