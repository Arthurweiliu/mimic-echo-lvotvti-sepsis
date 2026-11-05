## 96_fig1_shao_style.R — Figure 1 rebuilt in the BMC flowchart style of the published LVEF paper:
## rounded boxes, pale-blue cohort boxes with a steel-blue border, beige exclusion boxes connected by
## a horizontal spur, thick blue arrows, serif type, count in bold.
source("config.R")
suppressMessages(library(grid))

MAIN_FILL <- '#E8F0F7'; MAIN_BORDER <- '#2F5C88'
EXC_FILL  <- '#F2EDE3'; EXC_BORDER  <- '#B29674'
ARROW     <- '#3C6795'

png('Fig1_flow.png', width = 4.8, height = 7.4, units = 'in', res = 600)
grid.newpage()
pushViewport(viewport(width = 1, height = 1))

box <- function(x, y, w, h, lines, count, main = TRUE, cex_lab = 8.4, cex_cnt = 10) {
  grid.draw(roundrectGrob(x = unit(x, 'npc'), y = unit(y, 'npc'),
                          width = unit(w, 'npc'), height = unit(h, 'npc'),
                          r = unit(9, 'pt'),
                          gp = gpar(fill = if (main) MAIN_FILL else EXC_FILL,
                                    col = if (main) MAIN_BORDER else EXC_BORDER, lwd = 2.2)))
  n <- length(lines)
  # vertical layout: description lines above, count below
  span <- h * 0.72
  if (n == 1) {
    grid.text(lines[1], x = unit(x, 'npc'), y = unit(y + h * 0.20, 'npc'),
              gp = gpar(fontsize = cex_lab, fontfamily = 'serif', col = 'grey15'))
  } else {
    for (i in seq_along(lines)) {
      grid.text(lines[i], x = unit(x, 'npc'),
                y = unit(y + h * 0.30 - (i - 1) * (cex_lab / 72 / 7.4) * 1.35, 'npc'),
                gp = gpar(fontsize = cex_lab, fontfamily = 'serif', col = 'grey15'))
    }
  }
  grid.text(count, x = unit(x, 'npc'), y = unit(y - h * 0.22, 'npc'),
            gp = gpar(fontsize = cex_cnt, fontfamily = 'serif', fontface = 'bold', col = 'black'))
}

v_arrow <- function(x, y1, y2) {
  grid.lines(x = unit(c(x, x), 'npc'), y = unit(c(y1, y2), 'npc'),
             arrow = arrow(length = unit(7, 'pt'), type = 'closed'),
             gp = gpar(col = ARROW, lwd = 3.2, fill = ARROW))
}
spur <- function(x1, x2, y) {
  grid.lines(x = unit(c(x1, x2), 'npc'), y = unit(c(y, y), 'npc'),
             arrow = arrow(length = unit(6, 'pt'), type = 'closed'),
             gp = gpar(col = EXC_BORDER, lwd = 1.6, fill = EXC_BORDER))
}

mx <- 0.355; mw <- 0.615; mh <- 0.078
ex <- 0.815; ew <- 0.33;  eh <- 0.074
ys <- c(0.955, 0.826, 0.697, 0.568, 0.439, 0.310)

## main spine
box(mx, ys[1], mw, 0.092, c('Sepsis ICU admissions with a', 'transthoracic echocardiogram'),
    '14,136 admissions', cex_lab = 8.6, cex_cnt = 10.4)
box(mx, ys[2], mw, mh, c('One study-containing ICU stay', 'per patient'), '12,120 patients')
box(mx, ys[3], mw, mh, c('Structured LVOT-VTI', 'available'), '7,253 patients')
box(mx, ys[4], mw, mh, c('Index TTE 0-72 h after', 'ICU admission'), '5,366 patients')
box(mx, ys[5], mw, mh, c('Survived to the', 'index TTE'), '5,306 patients')
box(mx, ys[6], mw, mh, c('Analytic cohort: true first', 'ICU stay per patient'), '4,816 patients (1,573 deaths)',
    cex_cnt = 9.4)

## exclusions
mid <- function(i) (ys[i] + ys[i + 1]) / 2
box(ex, mid(1), ew, eh, c('Excluded - not the earliest', 'study-containing stay'), '2,016 admissions',
    main = FALSE, cex_lab = 7.4, cex_cnt = 8.4)
box(ex, mid(2), ew, eh, c('Excluded - no structured', 'LVOT-VTI'), '4,867 patients',
    main = FALSE, cex_lab = 7.4, cex_cnt = 8.4)
box(ex, mid(3), ew, eh, c('Excluded - index TTE', 'outside 0-72 h'), '1,887 patients',
    main = FALSE, cex_lab = 7.4, cex_cnt = 8.4)
box(ex, mid(4), ew, eh, c('Excluded - died before', 'the index TTE'), '60 patients',
    main = FALSE, cex_lab = 7.4, cex_cnt = 8.4)
box(ex, mid(5), ew, eh, c('Excluded - not the', 'first ICU stay'), '490 patients',
    main = FALSE, cex_lab = 7.4, cex_cnt = 8.4)

## connectors
for (i in 1:5) v_arrow(mx, ys[i] - mh / 2 - 0.004, ys[i + 1] + mh / 2 + 0.004)
for (i in 1:5) spur(mx + mw / 2 + 0.004, ex - ew / 2 - 0.006, mid(i))

## head-to-head subsets
h1x <- 0.225; h2x <- 0.625; hy <- 0.105; hw <- 0.37; hh <- 0.100
yb <- 0.232
grid.lines(x = unit(c(mx, mx), 'npc'), y = unit(c(ys[6] - mh / 2 - 0.004, yb), 'npc'),
           gp = gpar(col = ARROW, lwd = 3.0))
grid.lines(x = unit(c(h1x, h2x), 'npc'), y = unit(c(yb, yb), 'npc'), gp = gpar(col = ARROW, lwd = 3.0))
grid.lines(x = unit(c(h1x, h1x), 'npc'), y = unit(c(yb, hy + hh / 2 + 0.004), 'npc'),
           arrow = arrow(length = unit(7, 'pt'), type = 'closed'), gp = gpar(col = ARROW, lwd = 3.0, fill = ARROW))
grid.lines(x = unit(c(h2x, h2x), 'npc'), y = unit(c(yb, hy + hh / 2 + 0.004), 'npc'),
           arrow = arrow(length = unit(7, 'pt'), type = 'closed'), gp = gpar(col = ARROW, lwd = 3.0, fill = ARROW))
box(h1x, hy, hw, hh, c('H1: VTI + LVEF from', 'the same study'), '4,336 patients (1,403 deaths)',
    cex_lab = 7.8, cex_cnt = 8.8)
box(h2x, hy, hw, hh, c('H2: VTI + LVEF +', 'echo-derived CI'), '3,815 patients (1,220 deaths)',
    cex_lab = 7.8, cex_cnt = 8.8)

dev.off()
cat('Fig1 v8 saved\n')
