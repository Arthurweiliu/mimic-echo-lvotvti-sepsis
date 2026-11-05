## 99d_fig1_85mm_chain.R — Figure 1 at 85 mm, with the nested head-to-head subsets shown as a chain
## (4,816 -> H1 4,336 -> H2 3,815), not as a split.
source("config.R")
suppressMessages(library(grid))

W_IN <- 3.35; H_IN <- 5.90
BORDER <- 'black'; EXC_FILL <- '#EDEDED'

png('Fig1_flow.png', width = W_IN, height = H_IN, units = 'in', res = 600)
grid.newpage(); pushViewport(viewport(width = 1, height = 1))

geom <- function(n, cex_lab, cex_cnt, pad_in = 0.05) {
  line_h <- cex_lab * 1.26 / 72 / H_IN
  cnt_h  <- cex_cnt * 1.26 / 72 / H_IN
  gap    <- (line_h + cnt_h) / 2 + 0.020
  extent <- (n - 1) * line_h + gap + (line_h + cnt_h) / 2
  list(line_h = line_h, cnt_h = cnt_h, gap = gap, extent = extent,
       h = extent + 2 * pad_in / H_IN)
}
box <- function(x, yc, w, lines, count, excl = FALSE, cex_lab = 9, cex_cnt = 10) {
  n <- length(lines); g <- geom(n, cex_lab, cex_cnt)
  grid.draw(roundrectGrob(x = unit(x, 'npc'), y = unit(yc, 'npc'),
                          width = unit(w, 'npc'), height = unit(g$h, 'npc'),
                          r = unit(5, 'pt'),
                          gp = gpar(fill = if (excl) EXC_FILL else 'white',
                                    col = BORDER, lwd = if (excl) 0.9 else 1.5)))
  c0 <- yc - g$extent / 2 + g$cnt_h / 2
  for (i in seq_along(lines)) {
    grid.text(lines[i], x = unit(x, 'npc'),
              y = unit(c0 + g$gap + (n - i) * g$line_h, 'npc'),
              gp = gpar(fontsize = cex_lab, fontfamily = 'serif', col = 'grey10'))
  }
  grid.text(count, x = unit(x, 'npc'), y = unit(c0, 'npc'),
            gp = gpar(fontsize = cex_cnt, fontfamily = 'serif', fontface = 'bold'))
  g$h
}
arr <- function(x, y1, y2, lwd = 1.7) {
  grid.lines(x = unit(c(x, x), 'npc'), y = unit(c(y1, y2), 'npc'),
             arrow = arrow(length = unit(5, 'pt'), type = 'closed'),
             gp = gpar(col = 'black', lwd = lwd, fill = 'black'))
}
spur <- function(x1, x2, y) {
  grid.lines(x = unit(c(x1, x2), 'npc'), y = unit(c(y, y), 'npc'),
             arrow = arrow(length = unit(5, 'pt'), type = 'closed'),
             gp = gpar(col = 'black', lwd = 1.0, fill = 'black'))
}

mx <- 0.300; mw <- 0.545
ex <- 0.795; ew <- 0.385
y1 <- 0.915; y2 <- 0.735; y3 <- 0.555; y4 <- 0.360; y5 <- 0.165

h1 <- box(mx, y1, mw, c('Sepsis ICU admissions', 'with a TTE'), '14,136 admissions', cex_lab = 9.0, cex_cnt = 10.0)
h2 <- box(mx, y2, mw, c('First ICU stay with a', 'structured VTI'), '7,253 patients', cex_lab = 9.0, cex_cnt = 10.0)
h3 <- box(mx, y3, mw, c('Analytic cohort:', 'true first ICU stay'), '4,816 (1,573 deaths)', cex_lab = 9.0, cex_cnt = 9.6)

ex1y <- (y1 + y2) / 2; ex2y <- (y2 + y3) / 2
box(ex, ex1y, ew, c('No VTI 4,867;', 'not earliest (2,016)'), '6,883 excluded', excl = TRUE, cex_lab = 7.2, cex_cnt = 7.9)
box(ex, ex2y, ew, c('Outside 0-72 h 1,887;', 'died before 60;', 'not first stay 490'), '2,437 excluded', excl = TRUE, cex_lab = 7.2, cex_cnt = 7.9)

arr(mx, y1 - h1 / 2, y2 + h2 / 2)
arr(mx, y2 - h2 / 2, y3 + h3 / 2)
spur(mx, ex - ew / 2, ex1y)
spur(mx, ex - ew / 2, ex2y)

## nested head-to-head subsets: a chain, because H2 is a subset of H1
h4 <- box(mx, y4, mw, c('H1: VTI + LVEF from', 'the same study'), '4,336 (1,403 deaths)', cex_lab = 8.4, cex_cnt = 9.0)
h5 <- box(mx, y5, mw, c('H2: also a computable', 'echo-derived CI'), '3,815 (1,220 deaths)', cex_lab = 8.4, cex_cnt = 9.0)
arr(mx, y3 - h3 / 2, y4 + h4 / 2)
arr(mx, y4 - h4 / 2, y5 + h5 / 2)

## side notes for the two restrictions
note <- function(x, y, txt) {
  grid.text(txt, x = unit(x, 'npc'), y = unit(y, 'npc'), just = c('left', 'centre'),
            gp = gpar(fontsize = 7.2, fontfamily = 'serif', col = 'grey25', fontface = 'italic'))
}
note(0.335, (y3 + y4) / 2, '480 without same-study LVEF')
note(0.335, (y4 + y5) / 2, '521 without a computable eCI')

dev.off(); cat('Fig1 v12 (85 mm, nested chain) saved\n')
