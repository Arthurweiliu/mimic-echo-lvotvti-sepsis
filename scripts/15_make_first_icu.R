# 15_make_first_icu.R — 生成"仅首次 ICU"数据集；备份全量
setwd('F:/openclaw_workspace/sepsis_echo_lit')
d<-read.csv('h1_data.csv',stringsAsFactors=FALSE)
cc<-read.csv('cohort_covariates.csv',stringsAsFactors=FALSE)[,c('stay_id','icu_intime')]
m<-merge(d, cc, by='stay_id', all.x=TRUE)
m$t<-as.POSIXct(m$icu_intime)
m<-m[order(m$subject_id, m$t),]
first<-m[!duplicated(m$subject_id),]
# 备份全量，写入 first 为新的 h1_data.csv
if(!file.exists('h1_data_allstays.csv')) file.copy('h1_data.csv','h1_data_allstays.csv')
write.csv(first, 'h1_data_first.csv', row.names=FALSE)
cat('全量 stays=',nrow(d),' 患者=',length(unique(d$subject_id)),'\n')
cat('首次 ICU stays=',nrow(first),' 患者=',length(unique(first$subject_id)),' (每患者1条)\n')
cat('首次 ICU 中 VTI 可得=', sum(!is.na(as.numeric(first$vti))&as.numeric(first$vti)>0), '\n')
