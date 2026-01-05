#!/usr/bin/env Rscript

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(tools)


input_dir   <- "results/data"
pdf_panels  <- "results/data/all_smoothed_subplots.pdf"

files <- list.files(
  path = input_dir,
  pattern = "smoothed\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No *smoothed.csv files found in results/data/")
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

  # ---- Convert to long format ----
  df_fun <- tryCatch({
    bind_rows(
      tibble(x = df$x_loess, y = df$y_loess, method = "LOESS"),
      tibble(x = df$x_anwil, y = df$y_anwil, method = "ANWIL"),
      tibble(x = df$x_anwil_mode1, y = df$y_anwil_mode1, method = "ANWIL (mode 1)"),
      tibble(x = df$x_anwil_mode2, y = df$y_anwil_mode2, method = "ANWIL (mode 2)"),
      tibble(x = df$x_anwil_iterative, y = df$y_anwil_iterative, method = "ANWIL (iterative)"),
      tibble(x = df$x_nw, y = df$y_nw, method = "Nadaraya–Watson"),
      tibble(x = df$x_manle, y = df$y_manle, method = "ManLe"),
      tibble(x = df$x_amanle, y = df$y_amanle, method = "AmanLe")
    )
  }, error = function(e) {
    cat("Error processing long data:", e$message, "\n")
    NULL
  })

  if (is.null(df_fun)) next

  df_fun$method <- factor(df_fun$method, levels = c(
    "LOESS", "ANWIL", "ANWIL (mode 1)", "ANWIL (mode 2)", "ANWIL (iterative)", "Nadaraya–Watson", "ManLe", "AmanLe"
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
    geom_line(
      data = df_fun,
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

  print(p_fun)
}

dev.off()
cat("PDF subplots generated:", pdf_panels, "\n")