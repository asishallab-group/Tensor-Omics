# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(tools)

input_dir   <- "results/data"

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 9) {
  stop("Usage: Rscript plot_smooth_all.r <k_neighbors> <span> <max_iter> <kernel_type> <k_neighbors_sigma> <method_flag> <w_r> <w_e> <w_c>")
}

k_neighbors <- args[1]
span <- args[2]
max_iter <- args[3]
kernel_type <- args[4]
k_neighbors_sigma <- args[5]
method_flag <- args[6]
w_r <- sprintf("%.2f", as.numeric(args[7]))
w_e <- sprintf("%.2f", as.numeric(args[8]))
w_c <- sprintf("%.2f", as.numeric(args[9]))

# Print parameters to console
cat("k_neighbors:", k_neighbors, "\n")
cat("span:", span, "\n")
cat("max_iter:", max_iter, "\n")
cat("kernel_type:", kernel_type, "\n")
cat("k_neighbors_sigma:", k_neighbors_sigma, "\n")

# Define output path
pdf_panels <- paste0("results/plots/smoothed_k", k_neighbors, "_iter", max_iter, "_span", format(span, nsmall = 2), "_ksigma", k_neighbors_sigma, "_kernel", kernel_type, "_method", method_flag, "_wr", w_r, "_we", w_e, "_wc", w_c,".pdf")

# Update search pattern to include span
pattern <- paste0("smoothed_k", k_neighbors, "_iter", max_iter, "_span", span, "_ksigma", k_neighbors_sigma, "_kernel", kernel_type, "_method", method_flag, "_wr", w_r, "_we", w_e, "_wc", w_c, "\\.csv$")

files <- list.files(
  path = input_dir,
  pattern = pattern,
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No *smoothed.csv files found in results/data/")
}

# FIX 1: Rugosity (0..1) based on "following the general trend"
# ===============================================================
# Concept:
#   1) Build a robust trend of the original points (median by bins).
#   2) Compare each smooth curve against THIS trend (not raw points).
#   3) Measure "roughness" of the residual vs trend (abruptness + zig-zag).
#   4) Map to [0,1] using x/(x+k) (natural limit).
#
# Expected Result:
#   - ANWIL/wiggly curves => residual vs trend with high zig-zag => High Rugosity
#   - AmanLe/NW/LOESS (good results) => smooth residual vs trend => Low Rugosity

# ------------------------------------------------------------
# 1) Robust Trend (median by bins) from original points
# ------------------------------------------------------------
estimate_trend_binned <- function(x_ref, y_ref, nbins = 60) {
  ok <- complete.cases(x_ref, y_ref)
  x_ref <- x_ref[ok]; y_ref <- y_ref[ok]
  if (length(x_ref) < 10) return(NULL)

  brks <- seq(min(x_ref), max(x_ref), length.out = nbins + 1)
  bin  <- cut(x_ref, breaks = brks, include.lowest = TRUE)

  trend <- data.frame(x = x_ref, y = y_ref, bin = bin) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(
      x_med = median(x, na.rm = TRUE),
      y_med = median(y, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::filter(is.finite(x_med), is.finite(y_med)) |>
    dplyr::arrange(x_med)

  if (nrow(trend) < 5) return(NULL)
  trend
}

# ------------------------------------------------------------
# 2) Rugosity 0..1 on a series y(t) (scale-independent)
#    rough_amp = mean(|d2|)/mean(|d1|)  (relative abruptness)
#    zigzag    = frequency of sign changes in d2 (peaks/sharp shifts)
# ------------------------------------------------------------
calculate_rugosity01_series <- function(y, k = 0.25, w_amp = 0.7, w_zig = 0.3) {
  y <- y[is.finite(y)]
  if (length(y) < 6) return(NA_real_)

  d1 <- diff(y)
  d2 <- diff(d1)

  rough_amp <- mean(abs(d2), na.rm = TRUE) / (mean(abs(d1), na.rm = TRUE) + 1e-12)

  s2 <- sign(d2)
  s2 <- s2[s2 != 0]
  zigzag <- if (length(s2) < 2) 0 else sum(abs(diff(s2)) == 2) / length(s2)

  rough_raw <- w_amp * rough_amp + w_zig * zigzag

  # Stable mapping to [0,1]
  rough_raw / (rough_raw + k)
}

# ------------------------------------------------------------
# 3) Rugosity 0..1 of the RESIDUAL against the robust TREND
#    (This implements "Fix 1")
# ------------------------------------------------------------
calculate_rugosity01_residual_trend <- function(x_s, y_s, x_ref, y_ref,
                                                n = 200, nbins = 60,
                                                k = 0.25, w_amp = 0.7, w_zig = 0.3) {
  ok <- complete.cases(x_s, y_s)
  x_s <- x_s[ok]; y_s <- y_s[ok]
  if (length(x_s) < 5) return(NA_real_)

  trend <- estimate_trend_binned(x_ref, y_ref, nbins = nbins)
  if (is.null(trend)) return(NA_real_)

  # Uniform grid within smoother support
  xg <- seq(min(x_s), max(x_s), length.out = n)

  # Smooth curve on grid
  yhat <- approx(x_s, y_s, xout = xg, rule = 2)$y

  # Robust trend on the same grid
  ytr  <- approx(trend$x_med, trend$y_med, xout = xg, rule = 2)$y

  # Residual vs trend
  resid_trend <- ytr - yhat

  calculate_rugosity01_series(resid_trend, k = k, w_amp = w_amp, w_zig = w_zig)
}

# ------------------------------------------------------------
# 4) Bias 0..1 of residual vs trend (optional but useful)
#    Penalizes oversmoothing (systematic shift)
# ------------------------------------------------------------
calculate_bias01_residual_trend <- function(x_s, y_s, x_ref, y_ref,
                                            n = 200, nbins = 60, k = 0.05) {
  ok <- complete.cases(x_s, y_s)
  x_s <- x_s[ok]; y_s <- y_s[ok]
  if (length(x_s) < 5) return(NA_real_)

  trend <- estimate_trend_binned(x_ref, y_ref, nbins = nbins)
  if (is.null(trend)) return(NA_real_)

  xg <- seq(min(x_s), max(x_s), length.out = n)
  yhat <- approx(x_s, y_s, xout = xg, rule = 2)$y
  ytr  <- approx(trend$x_med, trend$y_med, xout = xg, rule = 2)$y
  resid_trend <- ytr - yhat

  sc <- IQR(ytr, na.rm = TRUE) + 1e-12
  b  <- abs(mean(resid_trend, na.rm = TRUE)) / sc

  b / (b + k)  # 0..1 mapping
}

# ------------------------------------------------------------
# 5) calculate_metrics using Fix 1 for Rugosity
# ------------------------------------------------------------
calculate_metrics <- function(df,
                              # Fix 1 parameters
                              trend_nbins = 5, grid_n = 200,
                              rug_k = 0.25, rug_w_amp = 0.3, rug_w_zig = 0.7,
                              bias_k = 0.05,
                              # score weights (heavy rugosity)
                              w_rmse = 0.35, w_rug = 0.50, w_cov = 0.10, w_bias = 0.05) {

  orig_x_base <- df$x_original
  orig_y_base <- df$y_original
  range_orig  <- diff(range(orig_x_base, na.rm = TRUE))

  name_map <- c(
    "loess" = "LOESS", "anwil" = "ANWIL",
    "anwil_iterative" = "ANWIL (iterative)", "nw" = "Nadaraya–Watson", "nw_knn" = "Nadaraya–Watson-Knn",
    "manle" = "ManLe", "amanle" = "AmanLe"
  )

  methods <- names(df)[grepl("^x_", names(df)) & names(df) != "x_original"]
  methods <- gsub("x_", "", methods)
  methods <- methods[methods != "manle_svd"]

  results <- lapply(methods, function(m) {
    x_s <- df[[paste0("x_", m)]]
    y_s <- df[[paste0("y_", m)]]
    ok <- !is.na(x_s) & !is.na(y_s)
    x_s <- x_s[ok]; y_s <- y_s[ok]
    if (length(x_s) < 5) return(NULL)

    # Maintain centering logic for ManLe/AmanLe
    if (m %in% c("manle", "amanle")) {
      target_x <- orig_x_base - mean(orig_x_base, na.rm = TRUE)
      target_y <- orig_y_base - mean(orig_y_base, na.rm = TRUE)
    } else {
      target_x <- orig_x_base
      target_y <- orig_y_base
    }

    # Coverage
    coverage <- diff(range(x_s)) / range_orig

    # RMSE in support intersection
    mask <- target_x >= min(x_s) & target_x <= max(x_s)
    y_interp <- approx(x_s, y_s, xout = target_x[mask], rule = 2)$y
    rmse <- sqrt(mean((target_y[mask] - y_interp)^2, na.rm = TRUE))

    # Fix 1: Rugosity (0..1) of residual vs robust trend
    rug01 <- calculate_rugosity01_residual_trend(
      x_s, y_s, target_x, target_y,
      n = grid_n, nbins = trend_nbins,
      k = rug_k, w_amp = rug_w_amp, w_zig = rug_w_zig
    )

    # (Optional) Bias vs trend (0..1)
    bias01 <- calculate_bias01_residual_trend(
      x_s, y_s, target_x, target_y,
      n = grid_n, nbins = trend_nbins, k = bias_k
    )

    data.frame(
      method     = ifelse(!is.na(name_map[m]), name_map[m], m),
      RMSE       = rmse,
      Coverage   = coverage,
      Rugosity   = rug01,
      Bias01     = bias01
    )
  })

  metrics_df <- bind_rows(results)

  # Normalization for RMSE and Coverage (0..1)
  norm01 <- function(v, reverse = FALSE) {
    mn <- min(v, na.rm = TRUE); mx <- max(v, na.rm = TRUE)
    if (!is.finite(mn) || !is.finite(mx) || mx == mn) return(rep(1, length(v)))
    res <- (v - mn) / (mx - mn)
    if (reverse) res <- 1 - res
    res
  }

  metrics_df$s_rmse <- norm01(metrics_df$RMSE, reverse = TRUE)
  metrics_df$s_cov  <- norm01(metrics_df$Coverage, reverse = FALSE)

  # Rugosity and Bias01 are already 0..1 (lower is better)
  metrics_df$s_rug  <- 1 - metrics_df$Rugosity
  metrics_df$s_bias <- 1 - metrics_df$Bias01

  metrics_df$Efficiency_Score <- round(
    (w_rmse * metrics_df$s_rmse +
     w_rug  * metrics_df$s_rug  +
     w_cov  * metrics_df$s_cov  +
     w_bias * metrics_df$s_bias) * 100, 1
  )

  metrics_df %>% arrange(desc(Efficiency_Score))
}

# ===========================================================
# 2) PDF SUBPLOTS (ONE PANEL PER METHOD)
# ===========================================================
pdf(pdf_panels, width = 12, height = 9)

for (f in files) {
  cat("Processing file:", f, "\n")

  df <- read_csv(f, show_col_types = FALSE)
  cat("Columns in the file:", colnames(df), "\n")

  if (!all(c("x_original", "y_original") %in% colnames(df))) {
    cat("Missing required columns in the file:", f, "\n")
    next
  }

  # ---- Convert to long format for plotting ----
  df_fun <- tryCatch({
    bind_rows(
      tibble(x = df$x_loess, y = df$y_loess, method = "LOESS"),
      tibble(x = df$x_anwil, y = df$y_anwil, method = "ANWIL"),
      tibble(x = df$x_anwil_iterative, y = df$y_anwil_iterative, method = "ANWIL (iterative)"),
      tibble(x = df$x_nw, y = df$y_nw, method = "Nadaraya–Watson"),
      tibble(x = df$x_nw, y = df$y_nw_knn, method = "Nadaraya–Watson-Knn"),
      tibble(x = df$x_manle, y = df$y_manle, method = "ManLe"),
      tibble(x = df$x_amanle, y = df$y_amanle, method = "AmanLe")
    )
  }, error = function(e) {
    cat("Error processing long data:", e$message, "\n")
    NULL
  })

  if (is.null(df_fun)) next

  df_fun$method <- factor(df_fun$method, levels = c(
    "LOESS", "ANWIL", "ANWIL (iterative)", "Nadaraya–Watson", "Nadaraya–Watson-Knn", "ManLe", "AmanLe"
  ))

  # DUPLICATE ORIGINAL POINTS FOR EACH METHOD FACET
  df_points_all <- df %>%
    select(x = x_original, y = y_original) %>%
    slice(rep(1:n(), times = length(unique(df_fun$method)))) %>%
    mutate(method = rep(unique(df_fun$method), each = nrow(df)))

  # Sort by x for correct line rendering
  df_fun <- df_fun %>% arrange(x)
  df_points_all <- df_points_all %>% arrange(x)

  # ManLe SVD points
  df_manle_extra <- df %>%
    select(x = x_manle_svd, y = y_manle_svd) %>%
    mutate(method = "ManLe")

  # Center values for ManLe/AmanLe facets
  df_points_all <- df_points_all %>%
    group_by(method) %>%
    mutate(
      x = ifelse(method %in% c("ManLe", "AmanLe"), x - mean(x[method %in% c("ManLe", "AmanLe")], na.rm = TRUE), x),
      y = ifelse(method %in% c("ManLe", "AmanLe"), y - mean(y[method %in% c("ManLe", "AmanLe")], na.rm = TRUE), y)
    ) %>%
    ungroup()

  # ---- Generate Visualization ----
  p_fun <- ggplot() +
    geom_point(
      data = df_points_all,
      aes(x = x, y = y),
      color = "black", alpha = 0.1, size = 1.1
    ) +
    geom_point(
      data = df_manle_extra,
      aes(x = x, y = y),
      color = "blue", alpha = 0.2, size = 1
    ) +
    geom_line(
      data = df_fun,
      aes(x = x, y = y, color = method),
      linewidth = 1
    ) +
    facet_wrap(~ method, scales = "free_y", ncol = 3) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste(file_path_sans_ext(basename(f))),
      subtitle = paste("k =", k_neighbors, "| span =", span, "| iter =", max_iter),
      x = "x",
      y = "y"
    ) 

  print(p_fun)
}

dev.off()
cat("PDF subplots generated:", pdf_panels, "\n")