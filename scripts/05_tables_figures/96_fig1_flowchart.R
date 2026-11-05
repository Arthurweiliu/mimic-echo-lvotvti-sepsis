## 98_fig1_final.R — Figure 1: three steps, two merged exclusion panels, full-width (170 mm),
## black and white, with visible connector lines (not just arrowheads) and >= 8 pt type.
source("config.R")
suppressMessages(library(grid))

W_IN <- 6.7     # 170 mm, the double-column width of the journal
H_IN <- 7.4
BORDER <- 'black'; EXC_FILL <- '#EDEDED'

png('Fig1_flow.png', width = W_IN, height = H_IN, units = 'in', res = 600)
grid.newpage(); pushViewport(viewport(width = 1, height = 1))

box <- function(x, y, w, h, lines, count, excl = FALSE, cex_lab = 9, cex_cnt = 10.5) {
  grid.draw(roundrectGrob(x = unit(x, 'npc'), y = unit(y, 'npc'),
                          width = unit(w, 'npc'), height = unit(h, 'npc'),
                          r = unit(10, 'pt'),
                          gp = gpar(fill = if (excl) EXC_FILL else 'white',
                                    col = BORDER, lwd = if (excl) 1.2 else 2.0)))
  n <- length(lines)
  step <- (cex_lab / 72) / H_IN * 1.55
  gap  <- 0.020                                   # space between last line and the count
  block <- (n - 1) * step + gap
  count_y <- y - block / 2 + 0.002
  for (i in seq_along(lines)) {
    grid.text(lines[i], x = unit(x, 'npc'), y = unit(count_y + gap + (n - i) * step, 'npc'),
              gp = gpar(fontsize = cex_lab, fontfamily = 'serif', col = 'grey10'))
  }
  grid.text(count, x = unit(x, 'npc'), y = unit(count_y, 'npc'),
            gp = gpar(fontsize = cex_cnt, fontfamily = 'serif', fontface = 'bold'))
}
arr <- function(x, y1, y2, lwd = 2.8) {
  grid.lines(x = unit(c(x, x), 'npc'), y = unit(c(y1, y2), 'npc'),
             arrow = arrow(length = unit(8, 'pt'), type = 'closed'),
             gp = gpar(col = 'black', lwd = lwd, fill = 'black'))
}
spur <- function(x1, x2, y) {   # from the vertical spine to the exclusion panel
  grid.lines(x = unit(c(x1, x2), 'npc'), y = unit(c(y, y), 'npc'),
             arrow = arrow(length = unit(8, 'pt'), type = 'closed'),
             gp = gpar(col = 'black', lwd = 1.5, fill = 'black'))
}

mx <- 0.315; mw <- 0.455
ex <- 0.780; ew <- 0.400
y1 <- 0.900; y2 <- 0.640; y3 <- 0.360
h1b <- 0.135; h2b <- 0.125; h3b <- 0.135; hx <- 0.135

box(mx, y1, mw, h1b, c('Sepsis ICU admissions with a', 'transthoracic echocardiogram'),
    '14,136 admissions', cex_lab = 9.5, cex_cnt = 11)
box(mx, y2, mw, h2b, c('First ICU stay with a', 'structured LVOT-VTI'), '7,253 patients',
    cex_lab = 9.5, cex_cnt = 11)
box(mx, y3, mw, h3b, c('Analytic cohort: true first', 'ICU stay per patient'),
    '4,816 patients (1,573 deaths)', cex_lab = 9.5, cex_cnt = 10.5)

ex1y <- (y1 + y2) / 2; ex2y <- (y2 + y3) / 2
box(ex, ex1y, ew, hx, c('Excluded: no structured LVOT-VTI', '(4,867), or not the earliest', 'study-containing stay (2,016)'),
    '6,883 admissions', excl = TRUE, cex_lab = 8.2, cex_cnt = 9.2)
box(ex, ex2y, ew, hx, c('Excluded: index study outside 0-72 h', '(1,887), died before it (60), or', 'not the first ICU stay (490)'),
    '2,437 patients', excl = TRUE, cex_lab = 8.2, cex_cnt = 9.2)

arr(mx, y1 - h1b / 2, y2 + h2b / 2)
arr(mx, y2 - h2b / 2, y3 + h3b / 2)
spur(mx, ex - ew / 2, ex1y)
spur(mx, ex - ew / 2, ex2y)

## head-to-head subsets
h1x <- 0.265; h2x <- 0.735; hy <- 0.120; hw <- 0.415; hh <- 0.095
yb <- 0.255
grid.lines(x = unit(c(mx, mx), 'npc'), y = unit(c(y3 - h3b / 2, yb), 'npc'), gp = gpar(col = 'black', lwd = 2.4))
grid.lines(x = unit(c(h1x, h2x), 'npc'), y = unit(c(yb, yb), 'npc'), gp = gpar(col = 'black', lwd = 2.4))
arr(h1x, yb, hy + hh / 2, lwd = 2.4)
arr(h2x, yb, hy + hh / 2, lwd = 2.4)
box(h1x, hy, hw, hh, c('H1: VTI + LVEF from', 'the same study'), '4,336 patients (1,403 deaths)',
    cex_lab = 8.6, cex_cnt = 9.2)
box(h2x, hy, hw, hh, c('H2: VTI + LVEF +', 'echo-derived CI'), '3,815 patients (1,220 deaths)',
    cex_lab = 8.6, cex_cnt = 9.2)

dev.off(); cat('Fig1 v10 saved\n')
