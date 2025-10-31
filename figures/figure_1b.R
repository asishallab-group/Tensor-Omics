require(circlize)
require(RColorBrewer)
require(latex2exp)

addAlpha <- function(col, alpha = .25) {
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha))
}

eucLen <- function(x) {
  sqrt(sum(sapply(x, function(i)
    i ^ 2)))
}

angleVecs <- function(x, y, conv.fac = (180 / pi)) {
  acos(sum(sapply(1:length(x), function(i)
    x[i] * y[i])) / (eucLen(x) * eucLen(y))) * conv.fac
}

colors <- brewer.pal(7, "Dark2")
alpha.cols <- sapply(colors, addAlpha)


dupl <- c(1.9, 1.2)
l.dupl <- eucLen(dupl)
orth <- c(.2, 1.75)
l.orth <- eucLen(orth)


pdf('figure_1b.pdf',
    width = 4,
    height = 4)

old.par <- par(mgp = c(.4, 0.4, 0))

plot(
  # rbind(orth, dupl),
  # pch = 18,
  # cex = 4,
  NA,
  xlim = c(0, 2.2),
  ylim = c(0, 2.2),
  xlab = TeX('\\textit{tissue} X'),
  # '',
  ylab = TeX('\\textit{tissue} Y'),
  # '',
  cex.lab = 1,
  bty = 'l',
  axes = FALSE
)
points(orth[1], orth[2], pch = 21, bg = alpha.cols[[1]], col = colors[[1]], cex = 2)
points(dupl[1], dupl[2], pch = 21, bg = alpha.cols[[2]], col = colors[[2]], cex = 2)



# lines(c(0,0,0,2.15), c(0,2.15,0,0))
draw.sector(
  45,
  45 + angleVecs(c(1, 1), orth),
  clock.wise = FALSE,
  center = c(0, 0),
  border = colors[[1]],
  col = alpha.cols[[1]],
  lwd = 2,
  rou1 = l.orth
)
draw.sector(
  45,
  45 - angleVecs(c(1, 1), dupl),
  clock.wise = TRUE,
  center = c(0, 0),
  border = colors[[2]],
  col = alpha.cols[[2]],
  lwd = 2,
  rou1 = l.dupl
)

# abline(0, 1, lwd = 4, col = colors[[4]])
# text(
#   .1,
#   2.15,
#   label = TeX('$\\textit{tissue} X$'),
#   pos = 3,
#   offset = .3,
#   cex = .75
# )
# text(
#   2.05,
#   0,
#   label = TeX('$\\textit{tissue} Y$'),
#   pos = 3,
#   offset = .5,
#   cex = .75
# )
text(
  orth[1],
  orth[2],
  label = TeX('$\\textit{o}$'),
  pos = 3,
  offset = 1,
  col = colors[[1]],
  cex = 1.5
)

text(
  orth[1],
  orth[2] + .17,
  label = TeX('$\\rightarrow$'),
  pos = 3,
  offset = 1,
  col = colors[[1]],
  cex = 1
)

text(
  orth[1],
  orth[2],
  label = TeX('$\\delta(\\textit{o})$'),
  pos = 4,
  offset = 2.25,
  col = colors[[1]],
  cex = 1.5
)

text(
  orth[1] + .19,
  orth[2] + .18,
  label = TeX('$\\rightarrow$'),
  pos = 4,
  offset = 2.25,
  col = colors[[1]],
  cex = 1
)

text(
  dupl[1],
  dupl[2],
  label = TeX('$\\textit{p}$'),
  pos = 1,
  offset = 1.25,
  col = colors[[2]],
  cex = 1.5
)

text(
  dupl[1],
  dupl[2] + .16,
  label = TeX('$\\rightarrow$'),
  pos = 1,
  offset = 1.25,
  col = colors[[2]],
  cex = 1
)

text(
  dupl[1] + .1,
  dupl[2],
  label = TeX('$\\delta(\\textit{p})$'),
  pos = 3,
  offset = 1.25,
  col = colors[[2]],
  cex = 1.5
)

text(
  dupl[1] + .17,
  dupl[2] + .22,
  label = TeX('$\\rightarrow$'),
  pos = 3,
  offset = 1.25,
  col = colors[[2]],
  cex = 1
)

text(
  1.85,
  1.65,
  label = TeX('$\\textit{D}$'),
  pos = 3,
  offset = 1.25,
  col = colors[[4]],
  cex = 1.5
)

text(
  1.85,
  1.87,
  label = TeX('$\\rightarrow$'),
  pos = 3,
  offset = 1.25,
  col = colors[[4]],
  cex = 1
)

# Calculate vectors
diff_vec <- dupl - orth
unit_vec <- diff_vec / eucLen(diff_vec)

# Difference vector
arrows(
  orth[1],
  orth[2],
  dupl[1],
  dupl[2],
  col = colors[[7]],
  lwd = 2.5,
  length = 0.15,
  cex = 1.5
)

# Unit vector
arrows(
  orth[1],
  orth[2] - .05,
  orth[1] + unit_vec[1],
  orth[2] + unit_vec[2] - .05,
  col = colors[[6]],
  lwd = 2.5,
  length = 0.1,
  cex = 1.5
)


# labels
text(
  (orth[1] + dupl[1]) / 2 + .2,
  (orth[2] + dupl[2]) / 2 - .1,
  label = TeX('$s_{\\textit{po}}$'),
  pos = 3,
  col = colors[[7]],
  cex = 1.5
)

text(
  (orth[1] + dupl[1]) / 2 + .14,
  (orth[2] + dupl[2]) / 2 + .14,
  label = TeX('\\rightarrow$'),
  pos = 3,
  col = colors[[7]],
  cex = 1
)

text(
  orth[1] + .4,
  orth[2] - .3,
  label = TeX('${s}_{\\textit{po}}^{norm}$'),
  pos = 1,
  col = colors[[6]],
  cex = 1.5
)

text(
  orth[1] + .2,
  orth[2] - .2,
  label = TeX('$\\rightarrow$'),
  pos = 1,
  col = colors[[6]],
  cex = 1
)

lines(c(0, 2), c(0, 2), col = colors[[4]], lwd = 2)
lines(c(0, 0), c(0, 2), col = "black", lwd = 2)
lines(c(0, 2), c(0, 0), col = "black", lwd = 2)

par(old.par)

dev.off()