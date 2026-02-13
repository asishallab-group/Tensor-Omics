#!/usr/bin/env Rscript

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(tools)


input_dir   <- "results/data"

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript plot_smooth_all.r <k_neighbors> <span> <n_iters_max> <kernel_type> <k_neighbors_sigma>")
}

k_neighbors <- as.integer(args[1])
span <- as.numeric(args[2])
n_iters_max <- as.integer(args[3])
kernel_type <- as.integer(args[4])
k_neighbors_sigma <- as.integer(args[5])

cat("k_neighbors:", k_neighbors, "\n")
cat("span:", span, "\n")
cat("n_iters_max:", n_iters_max, "\n")
cat("kernel_type:", kernel_type, "\n")
cat("k_neighbors_sigma:", k_neighbors_sigma, "\n")

pdf_panels <- paste0("results/plots/smoothed_k", k_neighbors, "_iter", n_iters_max, "_span", format(span, nsmall = 2), "_ksigma", k_neighbors_sigma, "_kernel", kernel_type, ".pdf")

# Update the pattern to include span
pattern <- paste0("smoothed_k", k_neighbors, "_iter", n_iters_max, "_span", format(span, nsmall = 2), "_ksigma", k_neighbors_sigma, "_kernel", kernel_type, "\\.csv$")

files <- list.files(
  path = input_dir,
  pattern = pattern,
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No *smoothed.csv files found in results/data/")
}

# FIX 1: Rugosidad (0..1) basada en "seguir la tendencia general"
# ===============================================================
# Idea:
#   1) Construimos una tendencia robusta de los puntos originales (mediana por bins).
#   2) Comparamos cada curva smooth contra ESA tendencia (no contra puntos crudos).
#   3) Medimos "aspereza" del residuo vs tendencia (brusquedad + zig-zag).
#   4) Mapeamos a [0,1] con x/(x+k) (límite natural).
#
# Resultado esperado:
#   - ANWIL/curvas wiggly => residuo vs tendencia con mucho zig-zag => Rugosity alta
#   - AmanLe/NW/LOESS buenos => residuo vs tendencia suave => Rugosity baja

library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# 1) Tendencia robusta (mediana por bins) a partir de puntos originales
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
# 2) Rugosidad 0..1 sobre una serie y(t) (sin depender de escala)
#    rough_amp = mean(|d2|)/mean(|d1|)  (brusquedad relativa)
#    zigzag    = frecuencia de cambios de signo en d2 (picos/cambios bruscos)
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

  # Mapeo estable a [0,1]
  rough_raw / (rough_raw + k)
}

# ------------------------------------------------------------
# 3) Rugosidad 0..1 del RESIDUO contra la TENDENCIA robusta
#    (Esto es el "fix 1")
# ------------------------------------------------------------
calculate_rugosity01_residual_trend <- function(x_s, y_s, x_ref, y_ref,
                                                n = 200, nbins = 60,
                                                k = 0.25, w_amp = 0.7, w_zig = 0.3) {
  ok <- complete.cases(x_s, y_s)
  x_s <- x_s[ok]; y_s <- y_s[ok]
  if (length(x_s) < 5) return(NA_real_)

  trend <- estimate_trend_binned(x_ref, y_ref, nbins = nbins)
  if (is.null(trend)) return(NA_real_)

  # Grid uniforme en el soporte del smoother
  xg <- seq(min(x_s), max(x_s), length.out = n)

  # Curva smooth en grid
  yhat <- approx(x_s, y_s, xout = xg, rule = 2)$y

  # Tendencia robusta en el mismo grid
  ytr  <- approx(trend$x_med, trend$y_med, xout = xg, rule = 2)$y

  # Residuo vs tendencia
  resid_trend <- ytr - yhat

  calculate_rugosity01_series(resid_trend, k = k, w_amp = w_amp, w_zig = w_zig)
}

# ------------------------------------------------------------
# 4) Bias 0..1 del residuo vs tendencia (opcional pero útil)
#    Castiga oversmoothing (desplazamiento sistemático)
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

  b / (b + k)  # 0..1
}

# ------------------------------------------------------------
# 5) calculate_metrics (tu función) usando Fix 1 para Rugosity
# ------------------------------------------------------------
calculate_metrics <- function(df,
                              # parámetros del fix 1
                              trend_nbins = 5, grid_n = 200,
                              rug_k = 0.25, rug_w_amp = 0.3, rug_w_zig = 0.7,
                              bias_k = 0.05,
                              # pesos del score (rugosidad pesada)
                              w_rmse = 0.35, w_rug = 0.50, w_cov = 0.10, w_bias = 0.05) {

  orig_x_base <- df$x_original
  orig_y_base <- df$y_original
  range_orig  <- diff(range(orig_x_base, na.rm = TRUE))

  name_map <- c(
    "loess" = "LOESS", "anwil" = "ANWIL",
    "anwil_mode1" = "ANWIL (mode 1)", "anwil_mode2" = "ANWIL (mode 2)",
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

    # Mantengo tu lógica de centrado para ManLe/AmanLe
    if (m %in% c("manle", "amanle")) {
      target_x <- orig_x_base - mean(orig_x_base, na.rm = TRUE)
      target_y <- orig_y_base - mean(orig_y_base, na.rm = TRUE)
    } else {
      target_x <- orig_x_base
      target_y <- orig_y_base
    }

    # Coverage
    coverage <- diff(range(x_s)) / range_orig

    # RMSE en intersección de soporte
    mask <- target_x >= min(x_s) & target_x <= max(x_s)
    y_interp <- approx(x_s, y_s, xout = target_x[mask], rule = 2)$y
    rmse <- sqrt(mean((target_y[mask] - y_interp)^2, na.rm = TRUE))

    # FIX 1: Rugosidad (0..1) del residuo vs tendencia robusta
    rug01 <- calculate_rugosity01_residual_trend(
      x_s, y_s, target_x, target_y,
      n = grid_n, nbins = trend_nbins,
      k = rug_k, w_amp = rug_w_amp, w_zig = rug_w_zig
    )

    # (Opcional) Bias vs tendencia (0..1)
    bias01 <- calculate_bias01_residual_trend(
      x_s, y_s, target_x, target_y,
      n = grid_n, nbins = trend_nbins, k = bias_k
    )

    data.frame(
      method     = ifelse(!is.na(name_map[m]), name_map[m], m),
      RMSE       = rmse,
      Coverage   = coverage,
      Rugosity = rug01,
      Bias01     = bias01
    )
  })

  metrics_df <- bind_rows(results)

  # Normalización simple para RMSE y Coverage (0..1) dentro del grupo
  norm01 <- function(v, reverse = FALSE) {
    mn <- min(v, na.rm = TRUE); mx <- max(v, na.rm = TRUE)
    if (!is.finite(mn) || !is.finite(mx) || mx == mn) return(rep(1, length(v)))
    res <- (v - mn) / (mx - mn)
    if (reverse) res <- 1 - res
    res
  }

  metrics_df$s_rmse <- norm01(metrics_df$RMSE, reverse = TRUE)
  metrics_df$s_cov  <- norm01(metrics_df$Coverage, reverse = FALSE)

  # Rugosity y Bias01 ya están en 0..1 (menor es mejor)
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

# ------------------------------------------------------------
# Tips rápidos de ajuste (sin tocar la lógica):
#   - Si la tendencia binned queda muy “dentada”, baja nbins (ej 30-40)
#   - Si no castiga suficiente wiggles, sube rug_w_zig (ej 0.5) y baja rug_w_amp
#   - Si Rugosity se te pega mucho en un valor, ajusta rug_k (ej 0.1..0.6)
# ------------------------------------------------------------




# ===========================================================
# 2) PDF SUBPLOTS (ONE PANEL PER METHOD)
# ===========================================================
pdf(pdf_panels, width = 12, height = 9)

# Evaluate metrics for each file
for (f in files) {
  cat("Processing file:", f, "\n")

  df <- read_csv(f, show_col_types = FALSE)
  cat("Columns in the file:", colnames(df), "\n")

  if (!all(c("x_original", "y_original") %in% colnames(df))) {
    cat("Missing required columns in the file:", f, "\n")
    next
  }

  # metrics <- calculate_metrics(df)
  # print(metrics)

  # # Save metrics to CSV
  # metrics_file <- paste0("results/metrics/metrics_", file_path_sans_ext(basename(f)), ".csv")
  # dir.create(dirname(metrics_file), showWarnings = FALSE, recursive = TRUE)
  # write_csv(metrics, metrics_file)
  # cat("Metrics saved to:", metrics_file, "\n")

  # ---- Convert to long format ----
  df_fun <- tryCatch({
    bind_rows(
      tibble(x = df$x_loess, y = df$y_loess, method = "LOESS"),
      tibble(x = df$x_anwil, y = df$y_anwil, method = "ANWIL"),
      tibble(x = df$x_anwil_mode1, y = df$y_anwil_mode1, method = "ANWIL (mode 1)"),
      tibble(x = df$x_anwil_mode2, y = df$y_anwil_mode2, method = "ANWIL (mode 2)"),
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
    "LOESS", "ANWIL", "ANWIL (mode 1)", "ANWIL (mode 2)", "ANWIL (iterative)", "Nadaraya–Watson", "Nadaraya–Watson-Knn", "ManLe", "AmanLe"
  ))

  # DUPLICATE ORIGINAL POINTS FOR EACH METHOD
  df_points_all <- df %>%
    select(x = x_original, y = y_original) %>%
    slice(rep(1:n(), times = length(unique(df_fun$method)))) %>%
    mutate(method = rep(unique(df_fun$method), each = nrow(df)))

  # Sort data by the x column before plotting
  df_fun <- df_fun %>% arrange(x)
  df_points_all <- df_points_all %>% arrange(x)

  # Add extra points for the "ManLe" method
  df_manle_extra <- df %>%
    select(x = x_manle_svd, y = y_manle_svd) %>%
    mutate(
      method = "ManLe"
    )

  # Center the x and y values for the points of the "ManLe" and "AmanLe" methods in df_points_all
  df_points_all <- df_points_all %>%
    group_by(method) %>%
    mutate(
      x = ifelse(method %in% c("ManLe", "AmanLe"), x - mean(x[method %in% c("ManLe", "AmanLe")], na.rm = TRUE), x),
      y = ifelse(method %in% c("ManLe", "AmanLe"), y - mean(y[method %in% c("ManLe", "AmanLe")], na.rm = TRUE), y)
    ) %>%
    ungroup()

  # ---- Generate the plot ----
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
    geom_point(
      data = df_fun %>% filter(grepl("^ANWIL", method)),
      aes(x = x, y = y, color = method),
      size = 0.5
    ) +
    geom_line(
      data = df_fun %>% filter(!grepl("^ANWIL", method)),
      aes(x = x, y = y, color = method),
      linewidth = 1
    ) +
    facet_wrap(~ method, scales = "free_y", ncol = 3) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste(file_path_sans_ext(basename(f))),
      x = "x",
      y = "y"
    ) 
    # +
    # geom_text(
    #   data = metrics,
    #   aes(x = -Inf, y = Inf, label = paste0("RMSE: ", RMSE, "\nMSE: ", round(RMSE^2, 4), "\nCoverage: ", Coverage, "\nRugosity: ", Rugosity, "\nScore: ", Efficiency_Score)),
    #   hjust = -0.1, vjust = 1.1, inherit.aes = FALSE,
    #   size = 3  # Reduce text size
    # )

  print(p_fun)
}

dev.off()
cat("PDF subplots generated:", pdf_panels, "\n")

