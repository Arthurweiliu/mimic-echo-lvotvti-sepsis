## 51_lactate_rebuild_v2.R — clean rebuild of the lactate supplement
## Primary covariate = first_day_bg.lactate_max (blood gas, [intime-6 h, intime+24 h]).
## Supplement source = clinical-chemistry lactate (labevents itemid 50813, serum),
##   max in [intime, intime+24 h]  ->  lac_lab24
## Fallback for sensitivity = nearest chemistry lactate within +/-12 h of that window.
source("config.R")
suppressMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(x))

d <- fread(file.path(DATA, "h1_data.csv"),
           colClasses=list(character=c('stay_id','subject_id')))
d[, `:=`(icu_intime = as.POSIXct(icu_intime, tz='UTC'), lactate = num(lactate))]
cat('h1_data rows', nrow(d), '| blood-gas lactate present', sum(!is.na(d$lactate)), '\n')

lab <- fread('lac.csv', select=c('subject_id','itemid','charttime','valuenum'))
lab <- lab[itemid==50813 & !is.na(valuenum)]
lab[, subject_id := as.character(subject_id)]
lab[, charttime := as.POSIXct(charttime, format='%Y/%m/%d %H:%M:%S', tz='UTC')]
lab <- lab[!is.na(charttime)]
cat('chemistry lactate rows', nrow(lab), '| subjects', uniqueN(lab$subject_id), '\n')

st <- d[, .(stay_id, subject_id=as.character(subject_id),
            w0=icu_intime, w1=icu_intime + 24*3600)]
setkey(st, subject_id, w0, w1)
lb <- copy(lab)[, `:=`(c0=charttime, c1=charttime)]; setkey(lb, subject_id, c0, c1)
ov <- foverlaps(lb, st, by.x=c('subject_id','c0','c1'), by.y=c('subject_id','w0','w1'), type='within')
win <- ov[!is.na(stay_id), .(lac_lab24=max(valuenum)), by=stay_id]

## nearest chemistry value within +/-12 h of the window
st2 <- st[, .(stay_id, subject_id, w0=w0 - 12*3600, w1=w1 + 12*3600)]
setkey(st2, subject_id, w0, w1)
ov2 <- foverlaps(lb, st2, by.x=c('subject_id','c0','c1'), by.y=c('subject_id','w0','w1'), type='within')
ov2 <- ov2[!is.na(stay_id)]
ov2 <- merge(ov2, d[, .(stay_id, icu_intime)], by='stay_id')
ov2[, dist_h := pmax(as.numeric(difftime(icu_intime, charttime, units='hours')),
                     as.numeric(difftime(charttime, icu_intime + 24*3600, units='hours')), 0)]
setorder(ov2, stay_id, dist_h, -valuenum)
near <- ov2[, .SD[1], by=stay_id][, .(stay_id, lac_near=valuenum, lac_near_dist_h=dist_h)]

d <- merge(d, win, by='stay_id', all.x=TRUE)
d <- merge(d, near, by='stay_id', all.x=TRUE)
d[, lactate_filled2 := fifelse(is.na(lactate), lac_lab24, lactate)]
d[, lactate_nearest := fifelse(is.na(lactate_filled2), lac_near, lactate_filled2)]
d[, lactate_src := fifelse(!is.na(lactate), 'blood gas',
                    fifelse(!is.na(lac_lab24), 'chemistry in window',
                     fifelse(!is.na(lac_near), 'chemistry nearest', 'missing')))]

cat('\n--- coverage (all 12,120) ---\n'); print(d[, .N, by=lactate_src])
cat('diabetes-check: chemistry in window but blood gas missing =',
    sum(is.na(d$lactate) & !is.na(d$lac_lab24)), '\n')

both <- d[!is.na(lactate) & !is.na(lac_lab24)]
cat('\nAgreement between blood-gas and chemistry lactate (both available, n =', nrow(both), '):\n')
cat('  Pearson r =', round(cor(both$lactate, both$lac_lab24), 4),
    '| Spearman r =', round(cor(both$lactate, both$lac_lab24, method='spearman'), 4),
    '| exact match =', sprintf('%.1f%%', 100*mean(abs(both$lactate-both$lac_lab24)<1e-9)),
    '| median diff =', round(median(both$lac_lab24-both$lactate), 3), '\n')
cat('  mean blood gas', round(mean(both$lactate),3), 'vs chemistry', round(mean(both$lac_lab24),3), '\n')

## nearest distances among patients rescued by the nearest rule
resc <- d[is.na(lactate_filled2) & !is.na(lac_near)]
cat('\nNearest-rule rescues:', nrow(resc), 'stays | median distance', round(median(resc$lac_near_dist_h),1),
    'h | IQR', paste(round(quantile(resc$lac_near_dist_h, c(.25,.75)),1), collapse='-'),
    '| values inside window', sum(resc$lac_near_dist_h==0), '\n')

fwrite(d[, .(stay_id, subject_id, lactate, lac_lab24, lac_near, lac_near_dist_h,
             lactate_filled2, lactate_nearest, lactate_src)],
       'lactate_v2.csv')
saveRDS(d[, .(stay_id, lactate, lac_lab24, lac_near, lac_near_dist_h,
              lactate_filled2, lactate_nearest, lactate_src)], 'lactate_v2.rds')
cat('\nwritten lactate_v2.csv\n')
