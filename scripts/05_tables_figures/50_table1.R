## 50_table1.R — Table 1 for the MI-primary cohort (n = 4,816) by VTI category
source("config.R")
suppressMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(x))

d <- fread(file.path(DATA, 'h1_data.csv'), colClasses=list(character=c('stay_id','subject_id')))
d[, icu_intime := as.POSIXct(icu_intime, tz='UTC')]
for (v in c('age','sofa','sapsii','charlson','lactate','creat','hgb','plt','wbc',
            'bili','inr','hr_mean','lvef','fu_days','event90','vti','vti_diam','norepi'))
  d[, (v) := num(get(v))]
d[, sex := ifelse(sex %in% c('M','F'), sex, NA)]
tf <- function(x) { s <- trimws(as.character(x)); s[is.na(s)] <- ''; ifelse(grepl('^t', tolower(s)), 1, 0) }
d[, `:=`(ventN = tf(vent), rrtN = tf(rrt))]

lk <- fread(file.path(DATA, "linked_cohort.csv"), colClasses=list(character='stay_id'))
lk <- lk[!duplicated(stay_id), .(stay_id, echo_dt_minus_intime_h)]
d <- merge(d, lk, by='stay_id', all.x=TRUE, sort=FALSE)
d[, echo_h := num(echo_dt_minus_intime_h)]
ic <- fread(file.path(DATA, "sepsis3_icustays.csv"), colClasses=list(character=c('subject_id','stay_id')))
ic[, intime := as.POSIXct(intime, tz='UTC')]
setorder(ic, subject_id, intime)
first_ids <- ic[!duplicated(subject_id), stay_id]
d[, start := echo_h/24]
A <- d[!is.na(vti) & vti>0 & !is.na(fu_days) & fu_days>0 & !is.na(echo_h) &
       echo_h>=0 & echo_h<=72 & fu_days>start & stay_id %in% first_ids]
A <- merge(A, fread('lactate_v2.csv', colClasses=list(character=c('stay_id','subject_id'))),
           by='stay_id', all.x=TRUE, suffixes=c('', '.v2'))
A[, lactate_primary := lactate_filled2]
A[, vcat := cut(vti, c(-Inf, 17, 20 + 1e-9, Inf), labels=c('<17', '17-20', '>20'), right=FALSE)]
cat('Table 1 cohort n =', nrow(A), '| deaths =', sum(A$event90), '\n')
print(table(A$vcat, useNA='ifany'))

cont <- list(
  list('Age, years', 'age'), list('SOFA', 'sofa'), list('SAPS-II', 'sapsii'),
  list('Charlson index', 'charlson'), list('Lactate, mmol/L', 'lactate_primary'),
  list('Creatinine, mg/dL', 'creat'), list('Haemoglobin, g/dL', 'hgb'),
  list('Platelets, K/uL', 'plt'), list('WBC, K/uL', 'wbc'),
  list('Bilirubin, mg/dL', 'bili'), list('INR', 'inr'),
  list('Heart rate, mean, bpm', 'hr_mean'), list('LVEF, %', 'lvef'),
  list('Time from ICU admission to index TTE, h', 'echo_h'))

fmt <- function(x) sprintf('%.1f (%.1f-%.1f)', median(x, na.rm=TRUE),
                           quantile(x, .25, na.rm=TRUE), quantile(x, .75, na.rm=TRUE))
groups <- c('<17', '17-20', '>20')
rows <- list()
for (r in cont) {
  x <- A[[r[[2]]]]
  cells <- c(fmt(x), sapply(groups, function(g) fmt(A[vcat==g][[r[[2]]]])))
  p <- tryCatch(kruskal.test(x ~ A$vcat)$p.value, error=function(e) NA)
  rows[[length(rows)+1]] <- data.table(Variable=r[[1]], Overall=cells[1], g1=cells[2], g2=cells[3],
                                       g3=cells[4], p=p, Missing=sum(is.na(x)))
}
bin <- list(
  list('Female, n (%)', function(z) z$sex=='F'),
  list('Vasopressor, n (%)', function(z) !is.na(z$norepi) & z$norepi>0),
  list('Invasive ventilation, n (%)', function(z) z$ventN==1),
  list('Renal replacement therapy, n (%)', function(z) z$rrtN==1),
  list('90-day mortality, n (%)', function(z) z$event90==1))
for (r in bin) {
  f <- r[[2]]; v <- f(A)
  pct <- function(vv) sprintf('%d (%.1f)', sum(vv, na.rm=TRUE), 100*mean(vv, na.rm=TRUE))
  cells <- c(pct(v), sapply(groups, function(g) pct(f(A[vcat==g]))))
  tb <- table(A$vcat, v)
  p <- tryCatch(suppressWarnings(chisq.test(tb)$p.value), error=function(e) NA)
  rows[[length(rows)+1]] <- data.table(Variable=r[[1]], Overall=cells[1], g1=cells[2], g2=cells[3],
                                       g3=cells[4], p=p, Missing=sum(is.na(v)))
}
T1 <- rbindlist(rows)
T1[, p := ifelse(is.na(p), 'NA', ifelse(as.numeric(p) < 0.001, '<0.001', formatC(as.numeric(p), format='f', digits=3)))]
setcolorder(T1, c('Variable','Overall','g1','g2','g3','p','Missing'))
n <- nrow(A)
OUT <- T1
setnames(OUT, c('Variable', sprintf('Overall (n = %s)', format(n, big.mark=',')),
                sprintf('< 17 cm (n = %s)', format(sum(A$vcat=='<17'), big.mark=',')),
                sprintf('17-20 cm (n = %s)', format(sum(A$vcat=='17-20'), big.mark=',')),
                sprintf('> 20 cm (n = %s)', format(sum(A$vcat=='>20'), big.mark=',')),
                'p', 'Missing, n'))

md <- c('### Table 1. Baseline characteristics by VTI category (primary MI cohort, n = 4,816)', '',
        paste0('| ', paste(names(OUT), collapse=' | '), ' |'),
        paste0('|', paste(rep('---', ncol(OUT)), collapse='|'), '|'))
for (i in 1:nrow(OUT)) md <- c(md, paste0('| ', paste(unlist(OUT[i]), collapse=' | '), ' |'))
md <- c(md, '',
  '*Continuous variables: median (IQR), Kruskal-Wallis test across the three VTI categories; categorical variables: n (%), chi-square test. Lactate is the first-day blood-gas value, supplemented where missing by the maximum serum chemistry value in the first 24 h (remaining missing values are listed and were handled by multiple imputation in the models; descriptive cells are not imputed). Ventilation and renal replacement therapy were recorded as absent where the source field was blank.*')
writeLines(md, 'Table1_mi_cohort.md')
fwrite(OUT, 'Table1_mi_cohort.csv')
cat('\n'); print(OUT)
