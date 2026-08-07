#!/usr/bin/env Rscript
# Reads the CSVs run_stc.py wrote for one dataset/parameter combination and produces a
# multi-page PDF report: clusters, their tangent spaces, and their intersections. The STC
# counterpart of branch `smoothing`'s r/plot_lomanle_spheres.R -- same visual language
# (points colored by region, circles for a local radius, segments for a tangent direction
# scaled by its singular value), but plotting what STC produces (ensembles, tangent bases,
# reconciled intersections) instead of what LoManLe produces (a stitched manifold skeleton).
#
# Deliberately depends on ggplot2 alone: this environment does not have ggforce, viridis, or
# patchwork installed, all of which the original script used (viridis only for its color
# scale, which ggplot2 has built in since 3.0 as scale_color_viridis_c/d; ggforce only for
# geom_circle, replaced below by a small polygon helper; patchwork only for optional
# side-by-side panels, dropped -- every plot is its own PDF page instead, as most of the
# original's own plots already were).
#
# Usage: Rscript plot_stc.R <prefix>
#   where <prefix> is what run_stc.py was given as --out-prefix, i.e. this script reads
#   <prefix>_{points,membership,ensembles,super_ensembles,super_ensembles_jsi,params}.{csv,json}
#   and writes <prefix>.pdf

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
membership_df <- read_or_empty(paste0(prefix, "_membership.csv"), point_id = integer(), ensemble_id = integer(), is_seed = integer())
ensembles_df <- read_or_empty(paste0(prefix, "_ensembles.csv"), ensemble_id = integer())
groups_df <- read_or_empty(paste0(prefix, "_super_ensembles.csv"), group_id = integer(), ensemble_id = integer())
jsi_df <- read_or_empty(paste0(prefix, "_super_ensembles_jsi.csv"),
                         group_id = integer(), ensemble_id_from = integer(), ensemble_id_to = integer(), jsi = numeric())

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
  keys <- c("k_min", "density_quantile", "alpha_max_deg", "d_max", "g_max", "f_max", "a", "o",
            "reconciliation_mode", "min_jsi", "n_vectors", "n_dimensions", "n_ensembles")
  keys <- keys[keys %in% names(params)]
  paste(sprintf("%s = %s", keys, unlist(params[keys])), collapse = ", ")
}
caption_text <- param_caption(params)

# Dimension columns are whatever is left in points_df after point_id/n_ensembles.
dim_names <- setdiff(names(points_df), c("point_id", "n_ensembles"))
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

# --- Page 7: Ensemble Reconciliation -- super-ensembles ---------------------------------
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
         subtitle = sprintf("mode=%s -- filled = grouped seeds, colored by super-ensemble; open circles = ungrouped ensembles",
                             if (!is.null(params)) params$reconciliation_mode else "?"),
         caption = caption_text)
  plot(p7)

  # --- Page 8: JSI along each super-ensemble's chain ---
  if (nrow(jsi_df) > 0) {
    from_xy <- merge(jsi_df, points_df[, c("point_id", x_name, y_name)],
                      by.x = "ensemble_id_from", by.y = "point_id", all.x = FALSE)
    # ensemble_id, not point_id -- re-merge properly via ensembles' own seed points
    from_xy <- merge(jsi_df, ensembles_df[, c("ensemble_id", "seed_point_id")],
                      by.x = "ensemble_id_from", by.y = "ensemble_id")
    from_xy <- merge(from_xy, points_df[, c("point_id", x_name, y_name)], by.x = "seed_point_id", by.y = "point_id")
    names(from_xy)[names(from_xy) == x_name] <- "x_from"
    names(from_xy)[names(from_xy) == y_name] <- "y_from"
    to_xy <- merge(jsi_df, ensembles_df[, c("ensemble_id", "seed_point_id")],
                    by.x = "ensemble_id_to", by.y = "ensemble_id")
    to_xy <- merge(to_xy, points_df[, c("point_id", x_name, y_name)], by.x = "seed_point_id", by.y = "point_id")
    names(to_xy)[names(to_xy) == x_name] <- "x_to"
    names(to_xy)[names(to_xy) == y_name] <- "y_to"
    edge_xy <- merge(from_xy[, c("group_id", "ensemble_id_from", "ensemble_id_to", "jsi", "x_from", "y_from")],
                      to_xy[, c("group_id", "ensemble_id_from", "ensemble_id_to", "x_to", "y_to")],
                      by = c("group_id", "ensemble_id_from", "ensemble_id_to"))

    p8 <- ggplot() +
      geom_point(data = points_df, aes(x = .data[[x_name]], y = .data[[y_name]]), color = "grey88", size = 0.4, alpha = 0.5) +
      geom_segment(data = edge_xy, aes(x = x_from, y = y_from, xend = x_to, yend = y_to, color = jsi),
                   linewidth = 1.1, alpha = 0.9) +
      geom_point(data = grouped_xy, aes(x = .data[[x_name]], y = .data[[y_name]]), color = "black", size = 1.2) +
      scale_color_viridis_c(option = "viridis", name = "JSI", limits = c(0, 1)) +
      coord_fixed() +
      theme_minimal() +
      labs(title = "Ensemble Reconciliation: Jaccard Similarity along each chain",
           subtitle = "One segment per consecutive pair within a super-ensemble's column (see mod_STC.md, super_ensembles_JSI)")
    plot(p8)
  }
} else {
  message("No super-ensembles found at this reconciliation threshold -- skipping the reconciliation pages.")
}

dev.off()
message(paste("Done:", pdf_out))
