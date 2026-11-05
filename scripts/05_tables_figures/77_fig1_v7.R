## 77_fig1_v7.R — flow figure for the MI cohort (H1 4,336; H2 3,815)
source("config.R")
png('Fig1_flow_v7.png', width=3600, height=3600, res=600)
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
box(0.305, 'H1: VTI + LVEF from the same study\n(n = 4,336; 1,403 deaths)', h=0.080)
box(0.165, 'H2: VTI + LVEF + echo-derived CI\n(n = 3,815; 1,220 deaths)', h=0.080)
arrows(mainx, ys[6]-0.030, mainx, 0.305+0.040, length=0.05, lwd=1.0)
arrows(mainx, 0.305-0.040, mainx, 0.165+0.040, length=0.05, lwd=1.0)
text(0.5, 1.005, 'Study flow', cex=1.05, font=2)
dev.off()
cat('Fig1 v7 saved\n')
