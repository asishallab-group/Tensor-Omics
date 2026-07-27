#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggforce)
  library(viridis)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript plot_lomanle_spheres.R <file.csv> <k_min> <g_threshold> <o_max> <o_min> [stability_threshold] [scale_factor] [max_iterations] [relative_conv_tol]")
}
# stability_threshold/scale_factor/max_iterations/relative_conv_tol are only
# used for figure captions, so older invocations without them still work.
for (idx in 6:9) {
  if (length(args) < idx) args[idx] <- NA
}

param_caption <- function(args, sep = ", ") {
  parts <- c(paste("k_min =", args[2]), paste("g_threshold =", args[3]),
             paste("o_max =", args[4]), paste("o_min =", args[5]))
  if (!is.na(args[6])) parts <- c(parts, paste("stability =", args[6]))
  if (!is.na(args[7])) parts <- c(parts, paste("scale_factor =", args[7]))
  if (!is.na(args[8])) parts <- c(parts, paste("max_iter =", args[8]))
  if (!is.na(args[9])) parts <- c(parts, paste("conv_tol =", args[9]))
  paste(parts, collapse = sep)
}

graphics.off()
while (dev.cur() > 1) dev.off()

df <- read.csv(args[1], check.names = FALSE)

# Skeleton backbone edges: lomanle now emits the actual line segments of the
# backbone (from its anchor-graph MST) in a companion "..._edges.csv" file, so
# no joining/ordering logic needs to happen here -- just read and draw them.
edges_file <- sub("\\.csv$", "_edges.csv", args[1])
if (file.exists(edges_file)) {
  edges_df <- read.csv(edges_file, check.names = FALSE)
  edges_df$stage <- ifelse(edges_df$stage == "iter1", "Iteration 1", "Final")
} else {
  message(paste("No companion edges file found at:", edges_file,
                 "-- backbone edges plot will be empty."))
  edges_df <- data.frame(stage = character(), edge_id = integer(),
                          x = numeric(), y = numeric(), xend = numeric(), yend = numeric())
}

# Normalize common coordinate column names if needed.
if (!("x" %in% names(df)) && ("x1" %in% names(df))) df$x <- df$x1
if (!("y" %in% names(df)) && ("x2" %in% names(df))) df$y <- df$x2
if (!("z" %in% names(df))) {
  if ("x3" %in% names(df)) df$z <- df$x3
  if (!("z" %in% names(df)) && ("y_original" %in% names(df)) && ("x1" %in% names(df)) && ("x2" %in% names(df))) {
    df$z <- df$y_original
  }
}

# Fill optional columns to keep plotting code robust.
if (!("n_anchors" %in% names(df))) df$n_anchors <- 1
if (!("label" %in% names(df))) df$label <- 0
if (!("anchor" %in% names(df))) df$anchor <- 0
if (!("radius" %in% names(df))) df$radius <- 0
if (!("density" %in% names(df))) df$density <- 0
if (!("gap" %in% names(df))) df$gap <- 0

if (!("sk_x" %in% names(df))) df$sk_x <- df$x
if (!("sk_y" %in% names(df))) df$sk_y <- df$y
if (!("sk_x_f" %in% names(df))) df$sk_x_f <- df$sk_x
if (!("sk_y_f" %in% names(df))) df$sk_y_f <- df$sk_y

if (("z" %in% names(df))) {
  if (!("sk_z" %in% names(df))) df$sk_z <- df$z
  if (!("sk_z_f" %in% names(df))) df$sk_z_f <- df$sk_z
}

if (!("v1_x" %in% names(df))) df$v1_x <- 0
if (!("v1_y" %in% names(df))) df$v1_y <- 0
if (!("v1_z" %in% names(df))) df$v1_z <- 0
if (!("s1" %in% names(df))) df$s1 <- 0

# Neighbor-selection diagnostics (see misc/smoothing_vecinos.md). Older CSVs
# without these columns fall back to values that make the new plots inert.
if (!("k_selected" %in% names(df))) df$k_selected <- 0
if (!("stability" %in% names(df))) df$stability <- 1
if (!("normal_error" %in% names(df))) df$normal_error <- 0
if (!("stopped_complex" %in% names(df))) df$stopped_complex <- 0

has_2d <- all(c("x", "y") %in% names(df))
has_3d <- all(c("x", "y", "z") %in% names(df))

if (!has_2d) {
  stop("Input CSV must contain at least 2D coordinates (x,y) or aliases (x1,x2).")
}

base_name <- gsub(".csv", "", basename(args[1]), fixed = TRUE)
out_dir <- "results/plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_2d_report <- function(df, edges_df, args, pdf_out) {
  message(paste("Saving 2D PDF report to:", pdf_out))
  pdf(pdf_out, width = 11, height = 8.5)

  message(sprintf("Edge rows loaded: iter1=%d final=%d",
                   sum(edges_df$stage == "Iteration 1"), sum(edges_df$stage == "Final")))
  if (nrow(edges_df) == 0) {
    message("No backbone edges found. The plot will show points only until the companion _edges.csv file is available.")
  }

  p1 <- ggplot(df) +
    geom_circle(aes(x0 = x, y0 = y, r = radius), color = "blue", alpha = 0.05, linewidth = 0.05) +
    geom_point(aes(x = x, y = y), size = 0.4, color = "black") +
    coord_fixed() +
    theme_minimal() +
    labs(title = "All Adaptive Spheres",
         subtitle = "Visualizing local neighborhood size (radius) for every point.")
  plot(p1)

  p2 <- ggplot(df) +
    geom_point(aes(x = x, y = y, color = density), size = 1.2) +
    scale_color_viridis(option = "magma") +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Point Density Distribution")
  plot(p2)

  p3 <- ggplot(df) +
    geom_point(aes(x = x, y = y), color = "grey90", size = 0.2, alpha = 0.5) +
    geom_segment(aes(x = x - v1_x * s1,
                     y = y - v1_y * s1,
                     xend = x + v1_x * s1,
                     yend = y + v1_y * s1,
                     color = density),
                 linewidth = 0.2, alpha = 0.4) +
    scale_color_viridis(option = "magma", name = "Density") +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Global Tangent Field (PC1)")
  plot(p3)

  anchors_df <- subset(df, anchor == 1)
  p4 <- ggplot(df) +
    geom_point(aes(x = x, y = y), color = "grey90", size = 0.3) +
    geom_circle(data = anchors_df,
                aes(x0 = x, y0 = y, r = radius, fill = density),
                color = "darkred", alpha = 0.2, linewidth = 0.3) +
    geom_point(data = anchors_df, aes(x = x, y = y), color = "black", size = 1) +
    scale_fill_viridis(option = "magma") +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Atlas Selection (Anchors)")
  plot(p4)

  anchor_svd_df <- subset(df, anchor == 1 & abs(v1_x) + abs(v1_y) > 0 & s1 > 0)
  p5 <- ggplot(df) +
    geom_point(aes(x = x, y = y), color = "grey80", size = 0.35, alpha = 0.6) +
    geom_segment(data = anchor_svd_df,
                 aes(x = x - v1_x * s1,
                     y = y - v1_y * s1,
                     xend = x + v1_x * s1,
                     yend = y + v1_y * s1),
                 color = "firebrick", linewidth = 0.55, alpha = 0.95) +
    geom_point(data = anchor_svd_df,
               aes(x = x, y = y),
               color = "firebrick", size = 1.6, alpha = 0.95) +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Anchor SVD Segments",
         subtitle = "Same anchors as plot 4: spheres first, then their SVD (PC1) in red")
  plot(p5)

  overlap_df <- subset(df, n_anchors >= 2 & label > 0)
  stitched_df <- overlap_df

  p6 <- ggplot(df) +
    geom_point(aes(x = x, y = y), color = "grey88", size = 0.4, alpha = 0.5) +
    geom_circle(data = anchors_df,
                aes(x0 = x, y0 = y, r = radius),
                color = "steelblue", fill = "steelblue", alpha = 0.08, linewidth = 0.3) +
    geom_point(data = overlap_df,
               aes(x = x, y = y, color = as.factor(label)),
               size = 2.0, alpha = 0.9) +
    geom_point(data = anchors_df,
               aes(x = x, y = y), color = "black", size = 1.5, shape = 17) +
    scale_color_viridis_d(option = "turbo", name = "Region") +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Intersection Regions (2+ anchors)",
    subtitle = "Anchor spheres in blue. Colored points = overlap zones. Triangles = anchor centers.",
         caption = paste("Parameters:", param_caption(args)))
  plot(p6)

  p7 <- ggplot(df) +
    geom_point(aes(x = x, y = y), color = "grey88", size = 0.4, alpha = 0.4) +
    geom_segment(data = stitched_df,
                 aes(x = x, y = y, xend = sk_x, yend = sk_y),
                 arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
                 color = "darkorange", alpha = 0.7, linewidth = 0.4) +
    geom_point(data = stitched_df,
               aes(x = sk_x, y = sk_y, color = as.factor(label)),
               size = 2.0, alpha = 0.95) +
    geom_point(data = anchors_df,
               aes(x = x, y = y), color = "black", size = 1.5, shape = 17) +
    scale_color_viridis_d(option = "turbo", name = "Region") +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Stitched Backbone (Iteration 1)",
         subtitle = "Orange path: stitching; Triangles: original anchors")
  plot(p7)

  single_df <- subset(df, n_anchors <= 1 | label == 0)
  junction_df <- subset(df, n_anchors >= 2 & label > 0)
  p8 <- ggplot() +
    geom_point(data = single_df,
               aes(x = sk_x_f, y = sk_y_f),
               color = "grey70", size = 1.5, alpha = 0.7) +
    geom_point(data = junction_df,
               aes(x = sk_x_f, y = sk_y_f, color = as.factor(label)),
               size = 1.5, alpha = 0.95) +
    geom_point(data = anchors_df,
               aes(x = x, y = y), color = "black", size = 1.5, shape = 17) +
    scale_color_viridis_d(option = "turbo", name = "Region") +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Final Stitched Backbone (Converged)",
         subtitle = "Final stitched points")
  plot(p8)

  p9 <- ggplot() +
    geom_point(data = df,
               aes(x = x, y = y),
               color = "grey75", size = 0.4, alpha = 0.45) +
    geom_segment(data = edges_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = "darkorange", linewidth = 1.0, alpha = 0.95,
                 lineend = "round") +
    facet_wrap(~stage, nrow = 1) +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Backbone Edges")
  plot(p9)

  # p11 <- ggplot(df) +
  #   geom_point(aes(x = x, y = y, color = stability), size = 1.2) +
  #   scale_color_viridis(option = "viridis", name = "Tangent\nstability") +
  #   coord_fixed() +
  #   theme_minimal() +
  #   labs(title = "Neighbor-Selection Diagnostics: Tangent Stability",
  #        subtitle = "Lower values mark points where growth stopped early because the tangent basis became unstable")
  # plot(p11)

  # complex_df <- subset(df, stopped_complex == 1)
  # p12 <- ggplot(df) +
  #   geom_point(aes(x = x, y = y), color = "grey85", size = 0.4, alpha = 0.5) +
  #   geom_point(data = complex_df, aes(x = x, y = y), color = "red", size = 1.4, alpha = 0.9) +
  #   coord_fixed() +
  #   theme_minimal() +
  #   labs(title = "Points Flagged as Geometrically Complex",
  #        subtitle = sprintf(
  #          "%d / %d points stopped growth early due to tangent instability (branches, intersections, curvature)",
  #          nrow(complex_df), nrow(df)))
  # plot(p12)

  # p13 <- ggplot(df) +
  #   geom_point(aes(x = x, y = y, color = k_selected), size = 1.2) +
  #   scale_color_viridis(option = "plasma", name = "k selected") +
  #   coord_fixed() +
  #   theme_minimal() +
  #   labs(title = "Neighbor-Selection Diagnostics: Selected Neighborhood Size")
  # plot(p13)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    moved_df <- df
    moved_df$disp <- sqrt((df$sk_x_f - df$sk_x)^2 + (df$sk_y_f - df$sk_y)^2)
    moved_df <- subset(moved_df, disp > 1e-6)

    p10a <- ggplot(df) +
      geom_point(aes(x = sk_x, y = sk_y,
                     color = as.factor(ifelse(n_anchors >= 2 & label > 0, label, 0))),
                 size = 0.8, alpha = 0.8) +
      scale_color_viridis_d(option = "turbo", guide = "none") +
      coord_fixed() + theme_minimal() +
      labs(title = "Iteration 1")

    p10b <- ggplot(df) +
      geom_point(aes(x = sk_x_f, y = sk_y_f,
                     color = as.factor(ifelse(n_anchors >= 2 & label > 0, label, 0))),
                 size = 0.8, alpha = 0.8) +
      scale_color_viridis_d(option = "turbo", guide = "none") +
      coord_fixed() + theme_minimal() +
      labs(title = "Converged")

    p10c <- ggplot(df) +
      geom_point(aes(x = sk_x, y = sk_y), color = "grey80", size = 0.5, alpha = 0.5) +
      geom_segment(data = moved_df,
                   aes(x = sk_x, y = sk_y, xend = sk_x_f, yend = sk_y_f),
                   arrow = arrow(length = unit(0.07, "cm"), type = "closed"),
                   color = "steelblue", alpha = 0.7, linewidth = 0.35) +
      coord_fixed() + theme_minimal() +
      labs(title = "Displacement iter1 -> final",
           subtitle = paste0("Points moved > 1e-6: ", nrow(moved_df), " / ", nrow(df)))

    p10 <- (p10a | p10b | p10c) +
      patchwork::plot_annotation(
        title = "Convergence: iteration 1 vs final",
        caption = paste("Parameters:", param_caption(args))
      )
    plot(p10)
  }

  dev.off()
}

plot_3d_report <- function(df, edges_df, args, html_out) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("3D plotting requires package 'plotly'. Install it with: install.packages('plotly')")
  }
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    stop("3D plotting requires package 'htmlwidgets'. Install it with: install.packages('htmlwidgets')")
  }

  message(paste("Saving 3D interactive report to:", html_out))

  anchors_df <- subset(df, anchor == 1)
  overlap_df <- subset(df, n_anchors >= 2 & label > 0)
  moved_df <- df
  moved_df$disp <- sqrt((df$sk_x_f - df$sk_x)^2 + (df$sk_y_f - df$sk_y)^2 + (df$sk_z_f - df$sk_z)^2)
  moved_df <- subset(moved_df, disp > 1e-6)

  # One consistent style per role, reused across every panel below: "Original"
  # (raw input, same data everywhere -- no more calling it "Background" in
  # some panels and "Original" in others) is a visible mid-grey, and every
  # other role gets its own solid, distinct color. Colors are fixed per role
  # (not mapped to a factor) so each shows up as exactly one legend entry --
  # mapping `color` to a factor in plotly splits it into one trace per level,
  # which is what produced duplicate-looking legend entries before.
  orig_marker  <- list(size = 3,   color = "rgba(90,90,90,0.55)")
  final_marker <- list(size = 3.5, color = "royalblue")
  anchor_marker <- list(size = 5, color = "black", symbol = "diamond")
  region_marker <- list(size = 4.5, color = "firebrick")

  p1 <- plotly::plot_ly(df, x = ~x, y = ~y, z = ~z,
                        type = "scatter3d", mode = "markers",
                        marker = orig_marker,
                        name = "Original") |>
    plotly::add_trace(data = df, x = ~sk_x_f, y = ~sk_y_f, z = ~sk_z_f,
                      type = "scatter3d", mode = "markers",
                      marker = final_marker, name = "Final stitched") |>
    plotly::add_trace(data = anchors_df, x = ~x, y = ~y, z = ~z,
                      type = "scatter3d", mode = "markers",
                      marker = anchor_marker,
                      name = "Anchors") |>
    plotly::layout(title = "3D: Original vs Final Stitched")

  p2 <- plotly::plot_ly(df, x = ~x, y = ~y, z = ~z,
                        type = "scatter3d", mode = "markers",
                        marker = orig_marker,
                        name = "Original") |>
    plotly::add_trace(data = overlap_df, x = ~x, y = ~y, z = ~z,
                      type = "scatter3d", mode = "markers",
                      marker = region_marker, name = "Intersection regions") |>
    plotly::layout(title = "3D: Intersection Regions")

  mk_segments <- function(d, x0, y0, z0, x1, y1, z1) {
    list(
      x = as.vector(rbind(d[[x0]], d[[x1]], NA_real_)),
      y = as.vector(rbind(d[[y0]], d[[y1]], NA_real_)),
      z = as.vector(rbind(d[[z0]], d[[z1]], NA_real_))
    )
  }

  p3 <- plotly::plot_ly(df, x = ~sk_x, y = ~sk_y, z = ~sk_z,
                        type = "scatter3d", mode = "markers",
                        marker = orig_marker,
                        name = "Iteration 1 points")
  if (nrow(moved_df) > 0) {
    seg_move <- mk_segments(moved_df, "sk_x", "sk_y", "sk_z", "sk_x_f", "sk_y_f", "sk_z_f")
    p3 <- p3 |>
      plotly::add_trace(x = seg_move$x, y = seg_move$y, z = seg_move$z,
                        type = "scatter3d", mode = "lines",
                        line = list(color = "steelblue", width = 3),
                        name = "Iter1 -> Final", showlegend = TRUE)
  }
  p3 <- p3 |>
    plotly::layout(title = "3D: Convergence Displacements")

  # Optional tangent vectors at anchors (if available)
  tangent_df <- subset(anchors_df, abs(v1_x) + abs(v1_y) + abs(v1_z) > 0 & s1 > 0)
  p4 <- plotly::plot_ly(df, x = ~x, y = ~y, z = ~z,
                        type = "scatter3d", mode = "markers",
                        marker = orig_marker,
                        name = "Original")
  if (nrow(tangent_df) > 0) {
    tangent_df$x0 <- tangent_df$x - tangent_df$v1_x * tangent_df$s1
    tangent_df$y0 <- tangent_df$y - tangent_df$v1_y * tangent_df$s1
    tangent_df$z0 <- tangent_df$z - tangent_df$v1_z * tangent_df$s1
    tangent_df$x1 <- tangent_df$x + tangent_df$v1_x * tangent_df$s1
    tangent_df$y1 <- tangent_df$y + tangent_df$v1_y * tangent_df$s1
    tangent_df$z1 <- tangent_df$z + tangent_df$v1_z * tangent_df$s1
    seg_tan <- mk_segments(tangent_df, "x0", "y0", "z0", "x1", "y1", "z1")
    p4 <- p4 |>
      plotly::add_trace(x = seg_tan$x, y = seg_tan$y, z = seg_tan$z,
                        type = "scatter3d", mode = "lines",
                        line = list(color = "firebrick", width = 4),
                        name = "PC1 segments", showlegend = TRUE)
  }
  p4 <- p4 |>
    plotly::layout(title = "3D: Anchor Tangent Field (PC1)")

  # Backbone skeleton in 3D -- same edges_df the 2D report's "Backbone Edges"
  # facet uses, just drawn as 3D line traces instead of ggplot segments.
  has_backbone_z <- all(c("z", "zend") %in% names(edges_df))
  edges_iter1_3d <- if (has_backbone_z) subset(edges_df, stage == "Iteration 1") else edges_df[0, ]
  edges_final_3d <- if (has_backbone_z) subset(edges_df, stage == "Final") else edges_df[0, ]

  p5 <- plotly::plot_ly(df, x = ~sk_x, y = ~sk_y, z = ~sk_z,
                        type = "scatter3d", mode = "markers",
                        marker = orig_marker,
                        name = "Iteration 1 points")
  if (nrow(edges_iter1_3d) > 0) {
    seg_iter1 <- mk_segments(edges_iter1_3d, "x", "y", "z", "xend", "yend", "zend")
    p5 <- p5 |>
      plotly::add_trace(x = seg_iter1$x, y = seg_iter1$y, z = seg_iter1$z,
                        type = "scatter3d", mode = "lines",
                        line = list(color = "darkorange", width = 5),
                        name = "Backbone", showlegend = TRUE)
  }
  p5 <- p5 |> plotly::layout(title = "3D: Backbone Skeleton (Iteration 1)")

  p6 <- plotly::plot_ly(df, x = ~sk_x_f, y = ~sk_y_f, z = ~sk_z_f,
                        type = "scatter3d", mode = "markers",
                        marker = final_marker,
                        name = "Final points")
  if (nrow(edges_final_3d) > 0) {
    seg_final <- mk_segments(edges_final_3d, "x", "y", "z", "xend", "yend", "zend")
    p6 <- p6 |>
      plotly::add_trace(x = seg_final$x, y = seg_final$y, z = seg_final$z,
                        type = "scatter3d", mode = "lines",
                        line = list(color = "darkorange", width = 5),
                        name = "Backbone", showlegend = TRUE)
  }
  p6 <- p6 |> plotly::layout(title = "3D: Backbone Skeleton (Final)")

  if (!has_backbone_z) {
    message("No z/zend columns in the edges file -- backbone skeleton panels will be empty.")
  }

  suppressWarnings({
    panel <- plotly::subplot(p1, p2, p3, p4, p5, p6, nrows = 2, margin = 0.04)
    panel <- panel |>
      plotly::layout(
        title = paste0("LOMANLE 3D report | ", param_caption(args, sep = " | "))
      )

    htmlwidgets::saveWidget(panel, file = html_out, selfcontained = TRUE)
  })
}

# Always generate 2D PDF projection for compatibility.
pdf_out <- file.path(out_dir, paste0(base_name, "_2d.pdf"))
plot_2d_report(df, edges_df, args, pdf_out)

# Generate 3D report when z exists.
if (has_3d) {
  html_out <- file.path(out_dir, paste0(base_name, "_3d.html"))
  plot_3d_report(df, edges_df, args, html_out)
  message("Done: generated both 2D and 3D outputs.")
} else {
  message("Done: generated 2D output. No z column detected, so 3D report was skipped.")
}
