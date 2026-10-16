# 35_hr_check.R — is HR itself adverse, and do SV (protective) and HR (adverse) cancel in CO?
{ .a <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
num <- function(x) suppressWarnings(as.numeric(x))
d <- read.csv('h1_data.csv', stringsAsFactors = FALSE)
for (v in c('age','sofa','charlson','lactate','creat','fu_days','event90','vti','vti_diam')) d[[v]] <- num(d[[v]])
d$sex <- factor(ifelse(d$sex %in% c('M','F'), d$sex, NA)); d$stay_id <- as.character(d$stay_id)
suppressMessages(library(survival))
cov <- 'age+sex+sofa+charlson+lactate+creat'
ad0 <- read.csv('analysis_data.csv', stringsAsFactors=FALSE, colClasses=c(stay_id='character'))
ad0 <- ad0[!duplicated(ad0$stay_id), c('stay_id','measurement_id','body_surface_area')]
d <- merge(d, ad0, by='stay_id', all.x=TRUE); d$bsa<-num(d$body_surface_area); d$measurement_id<-as.character(d$measurement_id)
ex <- read.csv('echo_extra.csv', stringsAsFactors=FALSE); ex$measurement_id<-as.character(ex$measurement_id)
hr <- ex[ex$measurement=='resting_hr', c('measurement_id','result')]; hr<-hr[!duplicated(hr$measurement_id),]; names(hr)[2]<-'rhr'
d <- merge(d, hr, by='measurement_id', all.x=TRUE); d$rhr<-num(d$rhr)
ccs <- c('event90','fu_days','age','sex','sofa','charlson','lactate','creat','vti')
A <- d[complete.cases(d[,ccs]) & num(d$fu_days)>0 & num(d$vti)>0, ]
A$sv <- A$vti*pi*(A$vti_diam/2)^2
A$hr <- ifelse(!is.na(A$rhr) & A$rhr>20 & A$rhr<220, A$rhr, NA)
H <- A[!is.na(A$sv) & A$sv>5 & A$sv<250 & !is.na(A$hr), ]
cat('N =', nrow(H), '\n')
g <- function(f,t,k=1){b<-summary(f)$coefficients; sprintf('HR=%.3f (%.3f-%.3f) p=%.2g',exp(b[t,1]*k),exp((b[t,1]-1.96*b[t,3])*k),exp((b[t,1]+1.96*b[t,3])*k),b[t,5])}
cat('HR  per +10 bpm (alone)      :', g(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ hr')),data=H),'hr',10), '\n')
cat('VTI per +5 cm  (alone)       :', g(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti')),data=H),'vti',5), '\n')
cat('SV  per +10 mL (alone)       :', g(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ sv')),data=H),'sv',10), '\n')
cat('VTI + HR (joint)  VTI        :', g(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti + hr')),data=H),'vti',5), '\n')
cat('VTI + HR (joint)  HR         :', g(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ vti + hr')),data=H),'hr',10), '\n')
cat('SV  + HR (joint)  SV         :', g(coxph(as.formula(paste('Surv(fu_days,event90)~',cov,'+ sv + hr')),data=H),'sv',10), '\n')
