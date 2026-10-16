# 16_flowchart.R — 研究流程图 (STROBE)，600 dpi
# Style: no background shading, black text only, per-line placement so text never
# touches the box borders.
# v2 (2026-10-14): adds the two downstream branches used by the revised analysis
#   Set B — co-primary landmark cohort (alive at 72 h)      n = 5,099
#   Set C — haemodynamic subset (SV, CO, CI computable)      n = 4,810
# --- portable working directory: 1st CLI argument, else this script's own directory ---
{ .a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  setwd(if (length(.a)) dirname(normalizePath(sub('^--file=', '', .a[1]))) else getwd()) }
png('H1_flowchart.png', 3600, 3900, res=600)
par(mar=c(0,0,0,0)); plot.new(); plot.window(c(0,15), c(-2.1,11.3))

CX <- 4.6          # centre of the main column
W  <- 8.6          # box width
LH <- 0.44         # line spacing (data units)

box <- function(y, lines, h, cex=0.95, font=1, w=W, cx=CX){
  rect(cx-w/2, y-h/2, cx+w/2, y+h/2, border='black', lwd=2, col='white')
  n <- length(lines)
  ys <- y + (n-1)/2*LH - (0:(n-1))*LH
  for(i in seq_along(lines)) text(cx, ys[i], lines[i], cex=cex, font=font)
}
down <- function(y1, y2){ arrows(CX, y1, CX, y2, length=0.12, lwd=2) }
excl <- function(y, lines){
  arrows(CX+W/2, y, 10.3, y, length=0.10, lwd=1.5)
  n <- length(lines); ys <- y + (n-1)/2*0.42 - (0:(n-1))*0.42
  for(i in seq_along(lines)) text(10.45, ys[i], lines[i], cex=0.88, adj=0)
}

y1 <- 10.40; y2 <- 9.05; y3 <- 7.50; y4 <- 5.85; y5 <- 4.35; y6 <- 2.95; y7 <- 1.45
h1 <- 1.00; h2 <- 1.00; h3 <- 1.38; h4 <- 1.38; h5 <- 1.00; h6 <- 1.00; h7 <- 1.25

box(y1, c("MIMIC-IV / MIMIC-IV-ECHO", "(adult ICU admissions)"), h1, font=2)
box(y2, c("Sepsis (Sepsis-3) ICU admissions", "n = 32,971"), h2)
box(y3, c("Transthoracic echocardiography", "within 72 h of ICU admission", "n = 14,136"), h3)
box(y4, c("One record per patient", "(first ICU stay retained)", "n = 12,120"), h4)
box(y5, c("LVOT-VTI available", "n = 7,253"), h5)
box(y6, c("Complete covariate data", "n = 5,410"), h6)
box(y7, c("Primary analysis cohort (Set A)", "n = 5,410 (1,823 deaths at 90 days)"), h7, font=2)

down(y1-h1/2, y2+h2/2); down(y2-h2/2, y3+h3/2); down(y3-h3/2, y4+h4/2)
down(y4-h4/2, y5+h5/2); down(y5-h5/2, y6+h6/2); down(y6-h6/2, y7+h7/2)

excl(y4, c("Subsequent ICU stays", "excluded (n = 2,016)"))
excl(y5, c("No LVOT-VTI", "(n = 4,867)"))
excl(y6, c("Missing covariates", "(n = 1,843)"))

## ---- downstream branches (Sets B and C) ----
y8 <- -0.85; h8 <- 1.75
LX <- 3.45; RX <- 10.1; WB <- 6.5
segments(CX, y7-h7/2, CX, y7-h7/2-0.28, lwd=2)
arrows(CX, y7-h7/2-0.28, LX, y8+h8/2, length=0.10, lwd=2)
arrows(CX, y7-h7/2-0.28, RX, y8+h8/2, length=0.10, lwd=2)
box(y8, c("Co-primary landmark cohort", "(Set B: alive at 72 h)", "n = 5,099 (1,512 deaths)",
          "311 died within 72 h"), h8, cex=0.82, w=WB, cx=LX)
box(y8, c("Haemodynamic subset", "(Set C: SV, CO, CI)", "n = 4,810",
          "600 unavailable"), h8, cex=0.82, w=WB, cx=RX)

dev.off(); cat('saved H1_flowchart.png\n')
