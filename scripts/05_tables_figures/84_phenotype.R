## 84_phenotype.R — Table S7: haemodynamic / geometric phenotype by VTI category
source("config.R")
suppressMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(x))

d <- fread(file.path(DATA, "h1_data.csv"), colClasses=list(character=c('stay_id','subject_id')))
for (v in c('vti','vti_diam','lvef','hr_mean','fu_days','event90','norepi','sofa','lactate'))
  d[, (v) := suppressWarnings(as.numeric(get(v)))]
lk <- fread(file.path(DATA, "linked_cohort.csv"), colClasses=list(character='stay_id'))
lk <- lk[!duplicated(stay_id), .(stay_id, echo_dt_minus_intime_h, measurement_id=as.character(measurement_id))]
d <- merge(d, lk, by='stay_id', all.x=TRUE, sort=FALSE)
d[, echo_h := num(echo_dt_minus_intime_h)]
ic <- fread(file.path(DATA, "sepsis3_icustays.csv"), colClasses=list(character=c('subject_id','stay_id')))
ic[, intime := as.POSIXct(intime, tz='UTC')]; setorder(ic, subject_id, intime)
first_ids <- ic[!duplicated(subject_id), stay_id]
d[, start := echo_h/24]
A <- d[!is.na(vti) & vti>0 & !is.na(fu_days) & fu_days>0 & !is.na(echo_h) & echo_h>=0 &
       echo_h<=72 & fu_days>start & stay_id %in% first_ids]

cov <- fread(file.path(DATA, "cohort_covariates.csv"), colClasses=list(character='stay_id'))
cov <- cov[!duplicated(stay_id), .(stay_id, mbp_mean=num(mbp_mean), uo=num(urineoutput_fd))]
A <- merge(A, cov, by='stay_id', all.x=TRUE, sort=FALSE)

ph <- fread(file.path(DATA, "tte_pheno.csv"),
            select=c('measurement_id','lvedd'))
ph[, measurement_id := as.character(measurement_id)]
ph <- ph[!duplicated(measurement_id)]
A <- merge(A, ph, by='measurement_id', all.x=TRUE, sort=FALSE)
A[, lvedd := num(lvedd)]

A[, sv := vti*pi*(vti_diam/2)^2]
A[, hr := ifelse(!is.na(hr_mean) & hr_mean > 25, hr_mean, NA)]
A[, co := sv*hr/1000]
A[, svr := ifelse(!is.na(mbp_mean) & !is.na(co) & co > 0, 80*mbp_mean/co, NA)]
A[, vcat := cut(vti, c(-Inf, 17, 20 + 1e-9, Inf), labels=c('<17','17-20','>20'), right=FALSE)]

vars <- list(
  list('Heart rate, bpm', 'hr', 1),
  list('Mean arterial pressure, mmHg', 'mbp_mean', 1),
  list('Stroke volume, mL', 'sv', 1),
  list('Cardiac output, L/min', 'co', 2),
  list('Estimated SVR, dyn s cm-5', 'svr', 0),
  list('LV end-diastolic dimension, cm', 'lvedd', 2),
  list('Norepinephrine, ug/kg/min', 'norepi', 3),
  list('Lactate, mmol/L', 'lactate', 1),
  list('Urine output, first 24 h, mL', 'uo', 0))

fmt <- function(x) sprintf('%.*f (%.*f-%.*f)', 0, median(x, na.rm=TRUE), 0,
                           quantile(x, .25, na.rm=TRUE), 0, quantile(x, .75, na.rm=TRUE))
fmt <- function(x, dg=1) sprintf(paste0('%.', dg, 'f (%.', dg, 'f-%.', dg, 'f)'),
                                 median(x, na.rm=TRUE), quantile(x, .25, na.rm=TRUE),
                                 quantile(x, .75, na.rm=TRUE))
rows <- list()
for (v in vars) {
  x <- A[[v[[2]]]]
  cells <- c(fmt(x, v[[3]]), sapply(c('<17','17-20','>20'), function(g) fmt(A[vcat==g][[v[[2]]]], v[[3]])))
  p <- tryCatch(kruskal.test(x ~ A$vcat)$p.value, error=function(e) NA)
  rows[[length(rows)+1]] <- data.frame(variable=v[[1]], overall=cells[1], g1=cells[2], g2=cells[3],
                                       g3=cells[4], p=p, missing=sum(is.na(x)))
}
T <- do.call(rbind, rows)
T$p <- ifelse(is.na(T$p), 'NA', ifelse(T$p < 0.001, '<0.001', formatC(T$p, format='f', digits=3)))
print(T, row.names=FALSE)
write.csv(T, 'TableS7_phenotype.csv', row.names=FALSE)

md <- c('### Table S7. Haemodynamic and geometric phenotype by VTI category', '',
        '| Variable | Overall | < 17 cm | 17-20 cm | > 20 cm | p | Missing, n |',
        '|---|---|---|---|---|---:|---:|')
for (i in seq_len(nrow(T))) md <- c(md, sprintf('| %s | %s | %s | %s | %s | %s | %d |',
  T$variable[i], T$overall[i], T$g1[i], T$g2[i], T$g3[i], T$p[i], T$missing[i]))
md <- c(md, '', paste0('*Median (IQR); Kruskal-Wallis test across the three categories. Cardiac output = stroke volume x echo heart rate / 1000; estimated SVR = 80 x mean arterial pressure / cardiac output. LV end-diastolic dimension was available only where the reporting template recorded it. Urine output is the first-24 h total from the derived covariates table.*'))
writeLines(md, 'TableS7_phenotype.md')
cat('\nwritten TableS7_phenotype.csv/.md\n')
