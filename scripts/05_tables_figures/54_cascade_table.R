## 54_cascade_table.R — assemble the missing-data cascade table from verified run output
source("config.R")
txt <- readLines('52_mi_primary_v2_out.txt', warn=FALSE)
pat <- '^(.+?)\\s+n=(\\d+) deaths=(\\d+) \\| HR ([0-9.]+) \\(([0-9.]+)-([0-9.]+)\\)(?: p=(\\S+))?(?: FMI=([0-9.]+))?'
m <- regmatches(txt, regexec(pat, txt))
rows <- Filter(function(z) length(z) >= 7, m)
if (!length(rows)) stop('no result lines parsed')
tab <- do.call(rbind, lapply(rows, function(z) {
  data.frame(analysis=trimws(z[2]), n=as.integer(z[3]), deaths=as.integer(z[4]),
             hr=as.numeric(z[5]), lo=as.numeric(z[6]), hi=as.numeric(z[7]),
             p=if (length(z) >= 8) z[8] else NA, fmi=if (length(z) >= 9) z[9] else NA,
             stringsAsFactors=FALSE)
}))
tab$estimand <- c('MI, m = 20', 'MI, m = 20', 'MI, m = 20', 'Complete case', 'Complete case', 'Complete case')[seq_len(nrow(tab))]
tab$lactate <- c('Blood gas only', 'Blood gas + chemistry (primary)', 'Blood gas + chemistry + nearest',
                 'Blood gas only', 'Blood gas + chemistry', 'Blood gas + chemistry + nearest')[seq_len(nrow(tab))]
tab$ci <- sprintf('%.2f (%.2f-%.2f)', tab$hr, tab$lo, tab$hi)
out <- tab[, c('estimand','lactate','n','deaths','ci','p','fmi')]
out$n <- format(as.integer(out$n), big.mark=',')
out$deaths <- format(as.integer(out$deaths), big.mark=',')
out$p[is.na(out$p) | out$p == ''] <- '-'
out$fmi[is.na(out$fmi) | out$fmi == ''] <- '-'
print(out)
write.csv(out, 'TableS_missing_data.csv', row.names=FALSE)
md <- c('### Table S. Missing-data sensitivity: VTI < 17 vs >= 17 cm for 90-day mortality',
        '',
        paste0('| ', paste(c('Analysis','Lactate source','n','Deaths','HR (95% CI)','p','FMI'), collapse=' | '), ' |'),
        '|---|---|---:|---:|---|---|---:|',
        apply(out, 1, function(r) paste0('| ', paste(r, collapse=' | '), ' |')),
        '',
        '*All models adjusted for age, sex, SOFA, Charlson index, lactate and creatinine, with left truncation at the index TTE. FMI, fraction of missing information. The primary analysis is multiple imputation (m = 20) in the 4,816-patient cohort; complete-case rows are shown for comparison.*')
writeLines(md, 'TableS_missing_data.md')
cat('\nwritten TableS_missing_data.md/.csv\n')
