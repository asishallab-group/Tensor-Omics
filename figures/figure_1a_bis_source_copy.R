require(circlize)
require(RColorBrewer)
require(latex2exp)

addAlpha <- function(col, alpha = .25) {
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha))
}

colors <- brewer.pal(7, "Dark2")
alpha.cols <- sapply(colors, addAlpha)

pdf('figure_1a_bis_source_copy.pdf', width = 8, height = 6)

plot(
  NA,
  pch = c(20, 18, 18),
  lwd = 2.5,
  cex = 3,
  xlim = c(0, 4),
  ylim = c(0, 6),  
  axes = FALSE,
  xlab = '',
  ylab = ''
)

for (i in 0:6){
    segments(x0 = i/2 + .6, y0 = .6, x1 = i/2 + .6, y1 = 4, col = "grey", lwd = .8)
}

for (i in 0:6){
    segments(x0 = .6, y0 = i/2 + .6, x1 = 4, y1 = i/2 + .6, col = "grey", lwd = .8)
}

lines(c(.6, .6), c(.6, 4), col = "black", lwd = 1.5)              
lines(c(.6, 4), c(.6, .6), col = "black", lwd = 1.5)              

# Tissue labels
text(0.2, 3.9, label = TeX('\\textit{leaf}'), pos = 4, offset = .2, cex = 1.0)
text(3.7, .3, label = TeX('\\textit{root}'), pos = 3, offset = .1, cex = 1.0)


centroid_1 <- c(1.5, 3.5)
centroid_2 <- c(3, 1.5)

points(centroid_1[1], centroid_1[2], pch = 21, bg = alpha.cols[[1]], col = colors[[1]], cex = 3)
points(centroid_2[1], centroid_2[2], pch = 21, bg = alpha.cols[[1]], col = colors[[1]], cex = 4)


centroid_1_pts <- list(
  p1 = centroid_1 + c(0.0, 0.4),
  p2 = centroid_1 + c(0.4, 0.1),
  p3 = centroid_1 + c(0.3, -0.15),
  p4 = centroid_1 + c(-0.2, -0.2),
  p5 = centroid_1 + c(-0.25, 0.2)
)

set.seed(1)

for (tp in names(centroid_1_pts)) {
  pt <- centroid_1_pts[[tp]]
  points(pt[1], pt[2], pch = 21, col = colors[sample(1:2, 1)], cex = 1.5)
}

centroid_2_pts <- list(
  p1 = centroid_2 + c(-0.25, -0.25),
  p2 = centroid_2 + c(0.2, 0.4),
  p3 = centroid_2 + c(0.3, -0.2),
  p4 = centroid_2 + c(-0.2, -0.5),
  p5 = centroid_2 + c(-0.32, 0.3)
)

for (tp in names(centroid_2_pts)) {
  pt <- centroid_2_pts[[tp]]
  points(pt[1], pt[2], pch = 21, col = colors[sample(1:2, 1)], cex = 2)
}

centroid_2_paras <- list(
    p1 = centroid_2 + c(-0.2, 1.2),
    p2 = centroid_2 + c(1.1, 0.3)
)

pt <- centroid_2_paras[[2]]
arrows(centroid_2_pts[[3]][1], centroid_2_pts[[3]][2], pt[1], pt[2],
        col = colors[2],
        lwd = 2.2,
        length = 0.1)
points(pt[1], pt[2], pch = 21, col = colors[2], cex = 2)

pt <- centroid_2_paras[[1]]
arrows(centroid_2_pts[[1]][1], centroid_2_pts[[1]][2], pt[1], pt[2],
        col = colors[2],
        lwd = 2.2,
        length = 0.1)
points(pt[1], pt[2], pch = 21, col = colors[2], cex = 2)

text(
  centroid_2_paras[[1]][1] - .25,
  centroid_2_paras[[1]][2] - .2,
  label = TeX('$\\textit{pea}$'),
  pos = 3,
  offset = 1,
  col = colors[[2]],
  cex = 1.2
)


text(
  centroid_2_paras[[2]][1] - .25,
  centroid_2_paras[[2]][2] - .15,
  label = TeX('$\\textit{bean}$'),
  pos = 3,
  offset = 1,
  col = colors[[2]],
  cex = 1.2
)


dev.off()