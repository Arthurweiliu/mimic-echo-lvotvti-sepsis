# 16_flowchart.R — 研究流程图 (STROBE)，600 dpi（加宽版）
setwd('F:/openclaw_workspace/sepsis_echo_lit')
png('H1_flowchart.png', 3600, 3400, res=600)
par(mar=c(0,0,0,0)); plot.new(); plot.window(c(0,15), c(0,11))
box <- function(y, txt, x=4.6, w=8.2, h=1.0, cex=0.9, font=1, col='white'){
  rect(x-w/2, y-h/2, x+w/2, y+h/2, border='black', lwd=2, col=col)
  text(x, y, txt, cex=cex, font=font)
}
down <- function(y1, y2){ arrows(4.6, y1, 4.6, y2, length=0.12, lwd=2) }
excl <- function(y, txt){ arrows(8.7, y, 10.6, y, length=0.1, lwd=1.5); text(10.7, y, txt, cex=0.82, adj=0) }

box(10.3, "MIMIC-IV / MIMIC-IV-ECHO\n(adult ICU admissions)", w=8.6, font=2)
box(8.9,  "Sepsis (Sepsis-3) ICU admissions\nn = 32,971")
box(7.5,  "Transthoracic echocardiography\nwithin 72 h of ICU admission\nn = 14,136")
box(6.1,  "One record per patient\n(first ICU stay retained)\nn = 12,120", col='#eef4ff')
box(4.7,  "LVOT-VTI available\nn = 7,253", col='#eef4ff')
box(3.2,  "Complete covariate data\nn = 5,410", col='#e7f7f4')
box(1.6,  "Primary analysis cohort\nn = 5,410 (1,823 deaths at 90 days)", w=8.8, h=1.15, font=2, col='#eaf7ec')
down(9.8, 9.4); down(8.4, 8.0); down(7.0, 6.6); down(5.6, 5.2); down(4.2, 3.7); down(2.7, 2.18)
excl(6.1, "Subsequent ICU stays\nexcluded (n = 2,016)")
excl(4.7, "No LVOT-VTI\n(n = 4,867)")
excl(3.2, "Missing covariates\n(n = 1,843)")
dev.off(); cat('saved H1_flowchart.png\n')
