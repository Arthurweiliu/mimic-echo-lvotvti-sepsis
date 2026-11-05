source("config.R")
suppressMessages(library(survival)); suppressMessages(library(splines))
z1 <- readRDS('mi_comp_H1.rds'); z2 <- readRDS('mi_comp_H2.rds')
X1 <- z1$X; X2 <- z2$X
cat('H1 n', nrow(X1), 'deaths', sum(X1$event90), '| mortality', sprintf('%.1f%%', 100*mean(X1$event90)), '\n')
cat('H2 n', nrow(X2), 'deaths', sum(X2$event90), '| mortality', sprintf('%.1f%%', 100*mean(X2$event90)), '\n')
cat('H1 VTI > 30:', sum(X1$vti > 30), '| > 32:', sum(X1$vti > 32), '\n')
cat('H1 VTI range', paste(range(X1$vti), collapse='-'), '| median', median(X1$vti), '\n')
cat('H2 eCI range', paste(round(range(X2$ci), 2), collapse='-'), '| median', round(median(X2$ci), 2), '\n')
cat('H1 LVEF median', median(X1$lvef), '\n')
## proportion of H1 in the >=25% averaged-risk category
p <- z1$avg[, 'clinical']
cat('H1 fraction with averaged risk >= 0.25:', sprintf('%.1f%%', 100*mean(p >= 0.25)), '\n')
