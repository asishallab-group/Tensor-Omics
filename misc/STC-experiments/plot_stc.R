#!/usr/bin/env Rscript
# Reads the CSVs run_stc.py wrote for one dataset/parameter combination and produces a
# multi-page PDF report (2D) and, for 3D+ data, an additional interactive HTML report:
# clusters, their tangent spaces, low-confidence fallback coverage, and which ensembles
# intersect. The STC counterpart of branch `smoothing`'s r/plot_lomanle_spheres.R -- same
# visual language (points colored by region, circles for a local radius, segments for a
# tangent direction scaled by its singular value), but plotting what STC produces (ensembles,
# tangent bases, low-confidence masks, reconciled intersections) instead of what LoManLe
# produces (a stitched manifold skeleton). Deliberately does NOT draw super-ensembles as an
# ordered "chain"/skeleton -- turning a set of intersecting ensembles into a single stitched
# manifold is LoManLe's job, not STC's; this script only ever reports which ensembles group
# together, never in what order or with what geometry they connect.
#
# Deliberately depends on ggplot2 alone for the 2D report: this environment does not have
# ggforce or patchwork installed, both of which the original script used (ggforce only for
# geom_circle, replaced below by a small polygon helper; patchwork only for optional
# side-by-side panels, dropped -- every plot is its own PDF page instead, as most of the
# original's own plots already were). viridis is not needed either way -- ggplot2 has had the
# same color scales built in since 3.0 as scale_color_viridis_c/d. The 3D report needs
# plotly + htmlwidgets (see plot_3d_report below); it is skipped, with a message, if either is
# unavailable.
#
# Usage: Rscript plot_stc.R <prefix>
#   where <prefix> is what run_stc.py was given as --out-prefix, i.e. this script reads
#   <prefix>_{points,membership,low_confidence_membership,ensembles,super_ensembles,
#   super_ensembles_overlap_coefficient,params}.{csv,json} and writes <prefix>.pdf, plus
#   <prefix>_3d.html when the data has 3 or more ambient dimensions.

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_stc.R <prefix>")
}
prefix <- args[1]

read_or_empty <- function(path, ...) {
  if (file.exists(path)) read.csv(path, check.names = FALSE) else data.frame(...)
}

points_df <- read.csv(paste0(prefix, "_points.csv"), check.names = FALSE)
if (!("n_low_confidence_ensembles" %in% names(points_df))) points_df$n_low_confidence_ensembles <- 0
membership_df <- read_or_empty(paste0(prefix, "_membership.csv"), point_id = integer(), ensemble_id = integer(), is_seed = integer())
low_confidence_df <- read_or_empty(paste0(prefix, "_low_confidence_membership.csv"), point_id = integer(), ensemble_id = integer())
overlap_coefficient_matrix_df <- read_or_empty(paste0(prefix, "_ensemble_overlap_coefficient_matrix.csv"), ensemble_id_1 = integer(), ensemble_id_2 = integer(), overlap_coefficient = numeric())
ensembles_df <- read_or_empty(paste0(prefix, "_ensembles.csv"), ensemble_id = integer())
groups_df <- read_or_empty(paste0(prefix, "_super_ensembles.csv"), group_id = integer(), ensemble_id = integer())

# run_stc.py's own params.txt: plain "key=value" lines, one per line -- read directly rather
# than via the sibling params.json, so this script needs no JSON-parsing package (jsonlite is
# not installed in every R environment this might run in, this one included).
read_params_txt <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path)
  lines <- lines[nzchar(lines)]
  parts <- strsplit(lines, "=", fixed = TRUE)
  values <- vapply(parts, function(p) paste(p[-1], collapse = "="), character(1))
  names(values) <- vapply(parts, `[`, character(1), 1)
  as.list(values)
}
params <- read_params_txt(paste0(prefix, "_params.txt"))

param_caption <- function(params) {
  if (is.null(params)) return("")
  keys <- c("k_min", "k_density", "bandwidth_percentile", "exclusion_radius_percentile",
            "chordal_dist_max_as_prcnt_of_range", "d_max", "g_max", "rmse_change_max", "f_max", "a", "o",
            "reconciliation_mode", "min_overlap_coefficient", "n_vectors", "n_dimensions", "n_ensembles")
  keys <- keys[keys %in% names(params)]
  paste(sprintf("%s = %s", keys, unlist(params[keys])), collapse = ", ")
}
caption_text <- param_caption(params)

# Dimension columns are whatever is left in points_df after point_id/n_ensembles.
dim_names <- setdiff(names(points_df), c("point_id", "n_ensembles", "n_low_confidence_ensembles"))
has_2d <- length(dim_names) >= 2
if (!has_2d) stop("Need at least 2 ambient dimensions to plot.")
x_name <- dim_names[1]
y_name <- dim_names[2]
has_3d <- length(dim_names) >= 3
z_name <- if (has_3d) dim_names[3] else NA

# A circle as a polygon (ggforce::geom_circle replacement, see the header note above).
circle_df <- function(cx, cy, r, id, n = 48) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(id = id, x = cx + r * cos(theta), y = cy + r * sin(theta))
}
circles_df <- function(cx, cy, r, ids, n = 48) {
  do.call(rbind, Map(function(a, b, c, d) circle_df(a, b, c, d, n), cx, cy, r, ids))
}

graphics.off()
while (dev.cur() > 1) dev.off()

pdf_out <- paste0(prefix, ".pdf")
message(paste("Saving PDF report to:", pdf_out))
pdf(pdf_out, width = 11, height = 8.5)

# --- Page 1: raw points, seeds highlighted ---------------------------------------------
seed_ids <- unique(membership_df$point_id[membership_df$is_seed == 1])
p1 <- ggplot(points_df, aes(x = .data[[x_name]], y = .data[[y_name]])) +
  geom_point(color = "grey70", size = 0.6, alpha = 0.7) +
  geom_point(data = subset(points_df, point_id %in% seed_ids), color = "black", size = 1.6, shape = 17) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Input point cloud and seeds",
       subtitle = sprintf("%d points, %d seed(s) (black triangles)", nrow(points_df), length(seed_ids)),
       caption = caption_text)
plot(p1)

# --- Page 2: clusters (ensemble membership), one color per ensemble --------------------
plot_df <- merge(points_df, membership_df, by = "point_id", all.x = FALSE)
unclustered_df <- subset(points_df, n_ensembles == 0)
p2 <- ggplot() +
  geom_point(data = unclustered_df, aes(x = .data[[x_name]], y = .data[[y_name]]),
             color = "grey85", size = 0.5, alpha = 0.6) +
  geom_point(data = plot_df, aes(x = .data[[x_name]], y = .data[[y_name]], color = as.factor(ensemble_id)),
             size = 1.0, alpha = 0.85) +
  scale_color_viridis_d(option = "turbo", name = "Ensemble", guide = if (nrow(ensembles_df) > 30) "none" else "legend") +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Clusters (ensembles)",
       subtitle = "Grey = never joined any ensemble. A point in 2+ ensembles is drawn once per ensemble, overplotted.",
       caption = caption_text)
plot(p2)

# --- Page 3: how many ensembles each point belongs to (overlap heat) -------------------
p3 <- ggplot(points_df, aes(x = .data[[x_name]], y = .data[[y_name]], color = n_ensembles)) +
  geom_point(size = 0.9) +
  scale_color_viridis_c(option = "magma", name = "# ensembles") +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Raw intersection density",
       subtitle = "How many ensembles each point belongs to, before Ensemble Reconciliation groups them")
plot(p3)

# --- Page 4: ensemble growth radii as spheres around each seed -------------------------
seeded_ens <- subset(ensembles_df, size > 0)
if (nrow(seeded_ens) > 0) {
  seed_xy <- merge(seeded_ens, points_df, by.x = "seed_point_id", by.y = "point_id")
  circ <- circles_df(seed_xy[[x_name]], seed_xy[[y_name]], seed_xy$growth_radius, seed_xy$ensemble_id)
  p4 <- ggplot() +
    geom_point(data = points_df, aes(x = .data[[x_name]], y = .data[[y_name]]), color = "grey85", size = 0.4) +
    geom_polygon(data = circ, aes(x = x, y = y, group = id), fill = "steelblue", color = "steelblue",
                 alpha = 0.06, linewidth = 0.2) +
    geom_point(data = seed_xy, aes(x = .data[[x_name]], y = .data[[y_name]]), color = "black", size = 1.2, shape = 17) +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Growth radii",
         subtitle = "One circle per ensemble, centered on its seed, radius = calc_ensemble_growth_radius")
  plot(p4)
}

# --- Page 5: tangent space, first principal direction per ensemble ---------------------
tangent_ens <- subset(ensembles_df, s1 > 0)
if (nrow(tangent_ens) > 0) {
  mu_x <- paste0("mu_", x_name)
  mu_y <- paste0("mu_", y_name)
  u1_x <- paste0("u1_", x_name)
  u1_y <- paste0("u1_", y_name)
  p5 <- ggplot() +
    geom_point(data = points_df, aes(x = .data[[x_name]], y = .data[[y_name]]), color = "grey85", size = 0.4, alpha = 0.6) +
    geom_segment(data = tangent_ens,
                 aes(x = .data[[mu_x]] - .data[[u1_x]] * s1, y = .data[[mu_y]] - .data[[u1_y]] * s1,
                     xend = .data[[mu_x]] + .data[[u1_x]] * s1, yend = .data[[mu_y]] + .data[[u1_y]] * s1,
                     color = as.factor(ensemble_id)),
                 linewidth = 0.8, alpha = 0.9) +
    geom_point(data = tangent_ens, aes(x = .data[[mu_x]], y = .data[[mu_y]], color = as.factor(ensemble_id)),
               size = 1.6) +
    scale_color_viridis_d(option = "turbo", guide = "none") +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Tangent space: first principal direction (U column 1)",
         subtitle = "One segment per ensemble, centered on mu, length scaled by tangent_scales' first entry")
  plot(p5)
}

# --- Page 6: second tangent direction, where d >= 2 -------------------------------------
tangent2_ens <- subset(ensembles_df, s2 > 0)
if (nrow(tangent2_ens) > 0) {
  mu_x <- paste0("mu_", x_name)
  mu_y <- paste0("mu_", y_name)
  u2_x <- paste0("u2_", x_name)
  u2_y <- paste0("u2_", y_name)
  p6 <- ggplot() +
    geom_point(data = points_df, aes(x = .data[[x_name]], y = .data[[y_name]]), color = "grey85", size = 0.4, alpha = 0.6) +
    geom_segment(data = tangent2_ens,
                 aes(x = .data[[mu_x]] - .data[[u2_x]] * s2, y = .data[[mu_y]] - .data[[u2_y]] * s2,
                     xend = .data[[mu_x]] + .data[[u2_x]] * s2, yend = .data[[mu_y]] + .data[[u2_y]] * s2,
                     color = as.factor(ensemble_id)),
                 linewidth = 0.8, alpha = 0.9) +
    scale_color_viridis_d(option = "turbo", guide = "none") +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Tangent space: second principal direction (U column 2)",
         subtitle = sprintf("%d / %d ensembles have an estimated intrinsic dimension >= 2", nrow(tangent2_ens), nrow(ensembles_df)))
  plot(p6)
} else {
  message("No ensemble has d >= 2 -- skipping the second-tangent-direction page.")
}

# --- Page 7: Ensemble Reconciliation -- super-ensembles (grouping only, no implied order) --
# Reports *which* ensembles intersect enough to group; deliberately does not draw them as a
# connected chain/skeleton -- see the header note on why that is out of scope here.
if (nrow(groups_df) > 0) {
  grouped_ens <- merge(groups_df, ensembles_df, by = "ensemble_id")
  grouped_xy <- merge(grouped_ens, points_df[, c("point_id", x_name, y_name)],
                       by.x = "seed_point_id", by.y = "point_id")
  ungrouped_ens <- subset(ensembles_df, size > 0 & !(ensemble_id %in% groups_df$ensemble_id))

  p7 <- ggplot() +
    geom_point(data = points_df, aes(x = .data[[x_name]], y = .data[[y_name]]), color = "grey88", size = 0.4, alpha = 0.5) +
    geom_point(data = ungrouped_ens_xy <- merge(ungrouped_ens, points_df[, c("point_id", x_name, y_name)],
                                                 by.x = "seed_point_id", by.y = "point_id"),
               aes(x = .data[[x_name]], y = .data[[y_name]]), color = "grey50", size = 1.4, shape = 1) +
    geom_point(data = grouped_xy, aes(x = .data[[x_name]], y = .data[[y_name]], color = as.factor(group_id)),
               size = 2.2, alpha = 0.9) +
    scale_color_viridis_d(option = "turbo", name = "Super-ensemble") +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Ensemble Reconciliation: super-ensembles",
         subtitle = sprintf("mode=%s -- filled = grouped seeds, colored by super-ensemble; open circles = ungrouped ensembles. Grouping only -- stitching a group into one ordered manifold is LoManLe's job, not STC's.",
                             if (!is.null(params)) params$reconciliation_mode else "?"),
         caption = caption_text)
  plot(p7)
} else {
  message("No super-ensembles found at this reconciliation threshold -- skipping the reconciliation page.")
}

# --- Page 7b: Overlap Coefficient heatmap, every non-empty ensemble pair (not just reconciled
# ones) -------------------------------------------------------------------------------------
if (nrow(overlap_coefficient_matrix_df) > 0) {
  ens_order <- sort(unique(c(overlap_coefficient_matrix_df$ensemble_id_1, overlap_coefficient_matrix_df$ensemble_id_2)))
  heat_df <- rbind(
    overlap_coefficient_matrix_df,
    data.frame(ensemble_id_1 = overlap_coefficient_matrix_df$ensemble_id_2, ensemble_id_2 = overlap_coefficient_matrix_df$ensemble_id_1,
               overlap_coefficient = overlap_coefficient_matrix_df$overlap_coefficient),
    data.frame(ensemble_id_1 = ens_order, ensemble_id_2 = ens_order, overlap_coefficient = 1.0)
  )
  p7b <- ggplot(heat_df, aes(x = as.factor(ensemble_id_1), y = as.factor(ensemble_id_2), fill = overlap_coefficient)) +
    geom_tile() +
    scale_fill_viridis_c(option = "viridis", name = "Overlap\nCoefficient", limits = c(0, 1)) +
    coord_fixed() +
    theme_minimal() +
    theme(axis.text = if (length(ens_order) > 40) element_blank() else element_text(size = 6),
          axis.ticks = element_blank()) +
    labs(title = "Ensemble Reconciliation: pairwise Overlap Coefficient heatmap",
         subtitle = sprintf("Every non-empty ensemble pair (%d ensembles) -- independent of reconciliation_mode/min_overlap_coefficient, unlike the super-ensembles page",
                             length(ens_order)),
         x = "ensemble_id", y = "ensemble_id")
  plot(p7b)
} else {
  message("Fewer than 2 non-empty ensembles -- skipping the Overlap Coefficient heatmap page.")
}

# --- Page 8: low-confidence fallback coverage -------------------------------------------
# Every point falls into exactly one of three coverage categories: in a *retained* ensemble
# (accepted, survived Stop Conditions); orphaned by the retained pipeline but covered by at
# least one seed's iteration-1 low-confidence mask (a fallback LoManLe may choose to use, see
# misc/mod_STC.md, "Ensemble identification", "Output"); or genuinely uncovered by anything.
coverage_category <- with(points_df, ifelse(
  n_ensembles > 0, "Retained ensemble",
  ifelse(n_low_confidence_ensembles > 0, "Low-confidence fallback only", "Uncovered")
))
points_df$coverage_category <- factor(coverage_category,
                                      levels = c("Retained ensemble", "Low-confidence fallback only", "Uncovered"))
coverage_counts <- table(points_df$coverage_category)
p8 <- ggplot(points_df, aes(x = .data[[x_name]], y = .data[[y_name]], color = coverage_category)) +
  geom_point(size = 0.9, alpha = 0.85) +
  scale_color_manual(values = c("Retained ensemble" = "steelblue", "Low-confidence fallback only" = "darkorange",
                                "Uncovered" = "grey70"),
                     name = "Coverage") +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Low-confidence fallback coverage",
       subtitle = sprintf("Retained=%d, low-confidence-only=%d, uncovered=%d -- orange points have no retained ensemble but a seed's iteration-1 mask still reaches them",
                           coverage_counts["Retained ensemble"], coverage_counts["Low-confidence fallback only"],
                           coverage_counts["Uncovered"]),
       caption = caption_text)
plot(p8)

# --- Page 9: parameters table (input this run used, plus estimate_stc_parameters' own
# proposals when requested) -- not a plot, a genuine table, via gridExtra::tableGrob. -------
if (!is.null(params) && length(params) > 0) {
  has_table_pkg <- requireNamespace("gridExtra", quietly = TRUE)
  if (!has_table_pkg) {
    message("Package 'gridExtra' is not installed -- skipping the parameters table page. Install it with ",
            "install.packages('gridExtra') (pure R, no system dependencies) to enable this.")
  } else {
    suppressPackageStartupMessages(library(gridExtra))

    shorten_path <- function(v) basename(as.character(v))
    all_keys <- names(params)
    input_keys <- all_keys[!grepl("^estimated_", all_keys)]
    estimated_keys <- all_keys[grepl("^estimated_", all_keys)]
    # Long absolute paths are noise in a printed table -- basename is enough to identify them.
    path_keys <- intersect(c("input_csv", "out_prefix"), input_keys)

    fmt_val <- function(key, value) {
      if (key %in% path_keys) return(shorten_path(value))
      num <- suppressWarnings(as.numeric(value))
      if (!is.na(num)) return(as.character(signif(num, 6)))
      value
    }

    # `keys` are the actual params[[...]] lookup keys; `labels` are what gets printed --
    # decoupled so the estimated table can strip its "estimated_" prefix for display while
    # still looking values up under their real (prefixed) key.
    build_table_df <- function(keys, labels = keys) {
      if (length(keys) == 0) return(NULL)
      data.frame(parameter = labels, value = vapply(keys, function(k) fmt_val(k, params[[k]]), character(1)),
                row.names = NULL)
    }
    input_table_df <- build_table_df(input_keys)
    estimated_table_df <- build_table_df(estimated_keys, sub("^estimated_", "", estimated_keys))

    table_theme <- gridExtra::ttheme_minimal(core = list(fg_params = list(hjust = 0, x = 0.05, cex = 0.8)),
                                             colhead = list(fg_params = list(fontface = "bold", cex = 0.85)))
    titled_table <- function(title, df) {
      gridExtra::arrangeGrob(grid::textGrob(title, gp = grid::gpar(fontface = "bold", cex = 0.9)),
                             gridExtra::tableGrob(df, rows = NULL, theme = table_theme),
                             heights = grid::unit(c(0.4, 9.6), "null"))
    }
    grobs <- list(titled_table("This run's own parameters", input_table_df))
    if (!is.null(estimated_table_df)) {
      grobs <- c(grobs, list(titled_table("estimate_stc_parameters (informational, never auto-applied)", estimated_table_df)))
    }
    grid::grid.newpage()
    gridExtra::grid.arrange(grobs = grobs, ncol = length(grobs), top = "Parameters",
                            bottom = "See misc/mod_STC.md, 'Estimate parameters from data' -- a heuristic starting point, not a converged answer")
  }
}

dev.off()
message(paste("Done:", pdf_out))

# --- 3D report, for 3+ ambient dimensions -------------------------------------------------
# The original intent here (see the smoothing branch's r/plot_lomanle_spheres.R,
# plot_3d_report) was an interactive plotly + htmlwidgets HTML export. That is not available
# in this environment: plotly/htmlwidgets/rgl all transitively need compiled system libraries
# (libcurl, libuv, openssl) whose *-devel headers are not installed here, installing them
# needs sudo, and this environment's own conda is broken (permission denied on the conda
# executable itself) -- all verified directly, not assumed. `scatterplot3d` is the fallback:
# pure R, no compiled system dependency beyond what R itself ships, so it installs and runs
# anywhere. It trades interactivity for a fixed viewing angle per page -- genuinely worse than
# the original plan, but a real fix for the actual bug (3D data silently only ever plotted in
# its first two dimensions), not a silent downgrade. If plotly/htmlwidgets ever become
# installable here (`sudo dnf install libcurl-devel openssl-devel libuv-devel`, then
# `install.packages(c("plotly","htmlwidgets"))`), this section is the one to replace.
if (has_3d) {
  has_scatterplot3d <- requireNamespace("scatterplot3d", quietly = TRUE)
  if (!has_scatterplot3d) {
    message("Package 'scatterplot3d' is not installed -- skipping the 3D report. Install it with ",
            "install.packages('scatterplot3d') (pure R, no system dependencies) to enable this.")
  } else {
    suppressPackageStartupMessages(library(scatterplot3d))

    pdf_3d_out <- paste0(prefix, "_3d.pdf")
    message(paste("Saving 3D PDF report to:", pdf_3d_out))
    pdf(pdf_3d_out, width = 9, height = 8)

    palette_for <- function(n) if (n < 1) character(0) else grDevices::hcl.colors(max(n, 3), "Viridis")[seq_len(n)]

    # --- 3D page 1: raw points, seeds highlighted ---
    s3d <- scatterplot3d(points_df[[x_name]], points_df[[y_name]], points_df[[z_name]],
                         color = "grey70", pch = 20, angle = 40,
                         xlab = x_name, ylab = y_name, zlab = z_name,
                         main = "Input point cloud and seeds (3D)",
                         sub = sprintf("%d points, %d seed(s) (black triangles). %s",
                                       nrow(points_df), length(seed_ids), caption_text))
    if (length(seed_ids) > 0) {
      seed_xyz <- subset(points_df, point_id %in% seed_ids)
      s3d$points3d(seed_xyz[[x_name]], seed_xyz[[y_name]], seed_xyz[[z_name]], col = "black", pch = 17)
    }

    # --- 3D page 2: clusters (ensemble membership), one color per ensemble ---
    if (nrow(ensembles_df) > 0) {
      ens_ids <- sort(unique(plot_df$ensemble_id))
      ens_colors <- setNames(palette_for(length(ens_ids)), ens_ids)
      s3d <- scatterplot3d(points_df[[x_name]], points_df[[y_name]], points_df[[z_name]],
                           color = "grey85", pch = 20, angle = 40,
                           xlab = x_name, ylab = y_name, zlab = z_name,
                           main = "Clusters (ensembles) (3D)",
                           sub = "Grey = never joined any ensemble. A point in 2+ ensembles is drawn once per ensemble, overplotted.")
      if (nrow(plot_df) > 0) {
        s3d$points3d(plot_df[[x_name]], plot_df[[y_name]], plot_df[[z_name]],
                     col = ens_colors[as.character(plot_df$ensemble_id)], pch = 19)
      }
    }

    # --- 3D page 3: low-confidence fallback coverage ---
    coverage_colors <- c("Retained ensemble" = "steelblue", "Low-confidence fallback only" = "darkorange",
                         "Uncovered" = "grey70")
    s3d <- scatterplot3d(points_df[[x_name]], points_df[[y_name]], points_df[[z_name]],
                         color = coverage_colors[as.character(points_df$coverage_category)], pch = 19, angle = 40,
                         xlab = x_name, ylab = y_name, zlab = z_name,
                         main = "Low-confidence fallback coverage (3D)",
                         sub = sprintf("Retained=%d, low-confidence-only=%d, uncovered=%d",
                                       coverage_counts["Retained ensemble"], coverage_counts["Low-confidence fallback only"],
                                       coverage_counts["Uncovered"]))
    legend("topright", legend = names(coverage_colors), col = coverage_colors, pch = 19, cex = 0.7, bty = "n")

    # --- 3D page 4: tangent space, first principal direction per ensemble ---
    tangent_ens <- subset(ensembles_df, s1 > 0)
    if (nrow(tangent_ens) > 0) {
      mu_x <- paste0("mu_", x_name); mu_y <- paste0("mu_", y_name); mu_z <- paste0("mu_", z_name)
      u1_x <- paste0("u1_", x_name); u1_y <- paste0("u1_", y_name); u1_z <- paste0("u1_", z_name)
      s3d <- scatterplot3d(points_df[[x_name]], points_df[[y_name]], points_df[[z_name]],
                           color = "grey85", pch = 20, angle = 40,
                           xlab = x_name, ylab = y_name, zlab = z_name,
                           main = "Tangent space: first principal direction (3D)",
                           sub = "One segment per ensemble, centered on mu, length scaled by tangent_scales' first entry")
      ens_colors <- setNames(palette_for(nrow(tangent_ens)), seq_len(nrow(tangent_ens)))
      for (i in seq_len(nrow(tangent_ens))) {
        row <- tangent_ens[i, ]
        p_from <- s3d$xyz.convert(row[[mu_x]] - row[[u1_x]] * row$s1, row[[mu_y]] - row[[u1_y]] * row$s1,
                                  row[[mu_z]] - row[[u1_z]] * row$s1)
        p_to <- s3d$xyz.convert(row[[mu_x]] + row[[u1_x]] * row$s1, row[[mu_y]] + row[[u1_y]] * row$s1,
                                row[[mu_z]] + row[[u1_z]] * row$s1)
        segments(p_from$x, p_from$y, p_to$x, p_to$y, col = ens_colors[i], lwd = 2)
      }
    }

    dev.off()
    message(paste("Done:", pdf_3d_out))
  }
} else {
  message("Data has fewer than 3 ambient dimensions -- no 3D report to generate.")
}
