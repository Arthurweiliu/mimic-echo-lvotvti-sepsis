# 15_make_first_icu.R — one record per patient (first ICU stay) → the analysis file h1_data.csv
# Pipeline: build_h1.mjs writes h1_data_allstays.csv (ALL eligible ICU stays); this script
# keeps the first ICU stay per patient and writes h1_data.csv (the file every downstream
# script reads) plus an identical h1_data_first.csv alias. No manual copying required.
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
d<-read.csv('h1_data_allstays.csv',stringsAsFactors=FALSE)
cc<-read.csv('cohort_covariates.csv',stringsAsFactors=FALSE)[,c('stay_id','icu_intime')]
m<-merge(d, cc, by='stay_id', all.x=TRUE)
m$t<-as.POSIXct(m$icu_intime)
m<-m[order(m$subject_id, m$t),]
first<-m[!duplicated(m$subject_id),]
# write the first-ICU dataset as h1_data.csv (the downstream analysis file) + alias
write.csv(first, 'h1_data_first.csv', row.names=FALSE)
write.csv(first, 'h1_data.csv', row.names=FALSE)
cat('全量 stays=',nrow(d),' 患者=',length(unique(d$subject_id)),'\n')
cat('首次 ICU stays=',nrow(first),' 患者=',length(unique(first$subject_id)),' (每患者1条)\n')
cat('首次 ICU 中 VTI 可得=', sum(!is.na(as.numeric(first$vti))&as.numeric(first$vti)>0), '\n')
