## 97_fig1_simple.R — Figure 1, black-and-white, three main steps and two merged exclusion panels
source("config.R")
suppressMessages(library(grid))

BORDER <- 'black'; EXC_FILL <- '#EDEDED'
LW_MAIN <- 1.8; LW_EXC <- 1.1

png('Fig1_flow.png', width = 4.7, height = 6.1, units = 'in', res = 600)
grid.newpage(); pushViewport(viewport(width = 1, height = 1))

box <- function(x, y, w, h, lines, count, excl = FALSE, cex_lab = 9, cex_cnt = 10.4) {
  grid.draw(roundrectGrob(x = unit(x, 'npc'), y = unit(y, 'npc'),
                          width = unit(w, 'npc'), height = unit(h, 'npc'),
                          r = unit(9, 'pt'),
                          gp = gpar(fill = if (excl) EXC_FILL else 'white',
                                    col = BORDER, lwd = if (excl) LW_EXC else LW_MAIN)))
  n <- length(lines); step <- (cex_lab / 72) / 6.0 * 1.45
  top <- if (n == 1) y + 0.012 else y + h * 0.30
  for (i in seq_along(lines)) {
    grid.text(lines[i], x = unit(x, 'npc'), y = unit(top - (i - 1) * step, 'npc'),
              gp = gpar(fontsize = cex_lab, fontfamily = 'serif', col = 'grey10'))
  }
  grid.text(count, x = unit(x, 'npc'), y = unit(y - h * 0.26, 'npc'),
            gp = gpar(fontsize = cex_cnt, fontfamily = 'serif', fontface = 'bold'))
}
arr <- function(x, y1, y2, lwd = 2.6) {
  grid.lines(x = unit(c(x, x), 'npc'), y = unit(c(y1, y2), 'npc'),
             arrow = arrow(length = unit(7, 'pt'), type = 'closed'),
             gp = gpar(col = 'black', lwd = lwd, fill = 'black'))
}
spur <- function(x1, x2, y) {
  grid.lines(x = unit(c(x1, x2), 'npc'), y = unit(c(y, y), 'npc'),
             arrow = arrow(length = unit(6, 'pt'), type = 'closed'),
             gp = gpar(col = 'black', lwd = 1.2, fill = 'black'))
}

mx <- 0.325; mw <- 0.48
ex <- 0.775; ew <- 0.385
y1 <- 0.905; y2 <- 0.640; y3 <- 0.360

box(mx, y1, mw, 0.140, c('Sepsis ICU admissions with a', 'transthoracic echocardiogram'),
    '14,136 admissions', cex_lab = 9.2, cex_cnt = 10.6)
box(mx, y2, mw, 0.130, c('First ICU stay with a', 'structured LVOT-VTI'), '7,253 patients',
    cex_lab = 9.2, cex_cnt = 10.6)
box(mx, y3, mw, 0.140, c('Analytic cohort: true first', 'ICU stay per patient'),
    '4,816 patients (1,573 deaths)', cex_lab = 9.2, cex_cnt = 10.2)

ex1y <- (y1 + y2) / 2; ex2y <- (y2 + y3) / 2
box(ex, ex1y, ew, 0.140, c('Excluded: no structured LVOT-VTI', '(4,867), or not the earliest', 'study-containing stay (2,016)'),
    '6,883 admissions', excl = TRUE, cex_lab = 7.4, cex_cnt = 8.6)
box(ex, ex2y, ew, 0.140, c('Excluded: study outside 0-72 h', '(1,887), died before it (60), or', 'not the first ICU stay (490)'),
    '2,437 patients', excl = TRUE, cex_lab = 7.4, cex_cnt = 8.6)

arr(mx, y1 - 0.070, y2 + 0.065)
arr(mx, y2 - 0.065, y3 + 0.070)
spur(mx + mw / 2 + 0.004, ex - ew / 2 - 0.008, ex1y)
spur(mx + mw / 2 + 0.004, ex - ew / 2 - 0.008, ex2y)

## head-to-head subsets
h1x <- 0.270; h2x <- 0.730; hy <- 0.115; hw <- 0.40; hh <- 0.095
yb <- 0.250
grid.lines(x = unit(c(mx, mx), 'npc'), y = unit(c(y3 - 0.070, yb), 'npc'), gp = gpar(col = 'black', lwd = 2.2))
grid.lines(x = unit(c(h1x, h2x), 'npc'), y = unit(c(yb, yb), 'npc'), gp = gpar(col = 'black', lwd = 2.2))
arr(h1x, yb, hy + hh / 2 + 0.004, lwd = 2.2)
arr(h2x, yb, hy + hh / 2 + 0.004, lwd = 2.2)
box(h1x, hy, hw, hh, c('H1: VTI + LVEF from', 'the same study'), '4,336 patients (1,403 deaths)',
    cex_lab = 8.0, cex_cnt = 8.6)
box(h2x, hy, hw, hh, c('H2: VTI + LVEF +', 'echo-derived CI'), '3,815 patients (1,220 deaths)',
    cex_lab = 8.0, cex_cnt = 8.6)

dev.off(); cat('Fig1 v9 saved\n')
