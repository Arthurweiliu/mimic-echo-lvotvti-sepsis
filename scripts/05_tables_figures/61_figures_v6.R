## 61_figures_v6.R — Fig 1 flow (MI cohort) and Fig 2 KM on the full eligible cohort
source("config.R")
suppressMessages(library(data.table)); suppressMessages(library(survival))
num <- function(x) suppressWarnings(as.numeric(x))

## ---- Fig 1 ------------------------------------------------------------
png('Fig1_flow_v6.png', width=3600, height=4400, res=600)
par(mar=c(0.5,0.5,0.5,0.5)); plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1.03))
mainx <- 0.36; w <- 0.54; sidex <- 0.76
box <- function(y, txt, h=0.065, cex=0.70) {
  rect(mainx-w/2, y-h/2, mainx+w/2, y+h/2, border='black', lwd=1.3)
  text(mainx, y, txt, cex=cex)
}
side <- function(y, txt) text(sidex, y, txt, cex=0.64, adj=0, col='grey20')
arr <- function(y1, y2) arrows(mainx, y1-0.030, mainx, y2+0.030, length=0.05, lwd=1.0)
ys <- c(0.925, 0.830, 0.735, 0.640, 0.545, 0.450)
txt <- c('Sepsis ICU admissions with\na transthoracic echocardiogram\n(n = 14,136)',
         'One ICU stay per patient\n(earliest TTE-containing stay)\n(n = 12,120)',
         'LVOT-VTI available\n(n = 7,253)',
         'Index TTE 0-72 h after ICU admission\n(n = 5,366)',
         'Survived to the index TTE\n(n = 5,306; 60 died before)',
         'True first ICU stay per patient\n(n = 4,816; 1,573 deaths)')
for (i in 1:6) { box(ys[i], txt[i], h=if (i == 1) 0.095 else 0.075); if (i < 6) arr(ys[i], ys[i+1]) }
for (i in 3:6) segments(mainx+w/2+0.01, ys[i], sidex-0.01, ys[i], lty=2, col='grey50')
side(ys[3], 'No LVOT-VTI:\n4,867')
side(ys[4], 'Outside 0-72 h:\n1,887')
side(ys[5], 'Died before TTE:\n60')
side(ys[6], 'Not first ICU stay:\n490')
box(0.315, 'H1: VTI + LVEF from the same study\n(n = 4,336; 1,403 deaths)', h=0.078)
box(0.195, 'H2: VTI + LVEF + echo-derived CI\n(n = 4,081; 1,311 deaths)', h=0.078)
box(0.065, 'Complete-covariate head-to-head subsets\nH1 n = 3,185 (1,113 deaths); H2 n = 2,781 (952 deaths)', h=0.082, cex=0.62)
arrows(mainx, ys[6]-0.030, mainx, 0.315+0.039, length=0.05, lwd=1.0)
arrows(mainx, 0.315-0.039, mainx, 0.195+0.039, length=0.05, lwd=1.0)
arrows(mainx, 0.195-0.039, mainx, 0.065+0.041, length=0.05, lwd=1.0)
text(0.5, 1.005, 'Study flow', cex=1.05, font=2)
dev.off()

## ---- Fig 2 ------------------------------------------------------------
d <- fread(file.path(DATA, "h1_data.csv"), colClasses=list(character=c('stay_id','subject_id')))
for (v in c('fu_days','event90','vti')) d[, (v) := num(get(v))]
lk <- fread(file.path(DATA, "linked_cohort.csv"), colClasses=list(character='stay_id'))
lk <- lk[!duplicated(stay_id), .(stay_id, echo_dt_minus_intime_h)]
d <- merge(d, lk, by='stay_id', all.x=TRUE, sort=FALSE)
d[, echo_h := num(echo_dt_minus_intime_h)]
ic <- fread(file.path(DATA, "sepsis3_icustays.csv"), colClasses=list(character=c('subject_id','stay_id')))
ic[, intime := as.POSIXct(intime, tz='UTC')]; setorder(ic, subject_id, intime)
first_ids <- ic[!duplicated(subject_id), stay_id]
d[, start := echo_h/24]
A <- as.data.frame(d[!is.na(vti) & vti>0 & !is.na(fu_days) & fu_days>0 & !is.na(echo_h) &
       echo_h>=0 & echo_h<=72 & fu_days>start & stay_id %in% first_ids])
A$grp <- cut(A$vti, c(-Inf,17,20+1e-9,Inf), labels=c('<17','17-20','>20'), right=FALSE)
A$time_ttE <- A$fu_days - A$start
cat('Fig2 cohort n', nrow(A), 'deaths', sum(A$event90), '\n'); print(table(A$grp))
sf <- survfit(Surv(time_ttE, event90) ~ grp, data=A)
png('Fig2_KM_VTI_v5.png', width=3000, height=2800, res=600)
par(mar=c(5.0,4.8,2.4,1.2), mgp=c(2.8,0.7,0))
plot(sf, col=c('#b2182b','#2166ac','#4d9221'), lwd=2.6, xlim=c(0,90), ylim=c(0.4,1.0),
     xlab='Days from the index TTE', ylab='Survival probability', bty='n', las=1, mark.time=FALSE)
legend('bottomleft', c('VTI < 17 cm','VTI 17-20 cm','VTI > 20 cm'),
       col=c('#b2182b','#2166ac','#4d9221'), lwd=2.6, bty='n', cex=0.9)
dev.off()
cat('figures v6 saved\n')
