# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(tools)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop("Usage: Rscript plot_selected_methods.R <k_neighbors> <span> <max_iter> <kernel_type> <k_neighbors_sigma> <method_flag> <w_r> <w_e> <w_c>")
}

k <- args[1]
span <- args[2]
max_iter <- args[3]
kernel_type <- args[4]
k_neighbors_sigma <- args[5]
method_flag <- args[6]
w_r <- sprintf("%.2f", as.numeric(args[7]))
w_e <- sprintf("%.2f", as.numeric(args[8]))
w_c <- sprintf("%.2f", as.numeric(args[9]))

input_dir <- "results/data/2d"
output_dir <- "results/plots"
output_pdf <- file.path(output_dir, paste0("combined_methods_k", k, "_iter", max_iter, "_span", span, "_ksigma", k_neighbors_sigma, "_kernel", kernel_type, "_method", method_flag, "_wr", w_r, "_we", w_e, "_wc", w_c, "_nw.pdf"))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Define search patterns for main, sigma, and score history files
pattern_main <- paste0("smoothed_k", k, "_iter", max_iter, "_span", span, "_ksigma", k_neighbors_sigma, "_kernel", kernel_type, "_method", method_flag, "_wr", w_r, "_we", w_e, "_wc", w_c, "\\.csv$")
pattern_sigma <- "_anwil_std\\.csv$"
pattern_score_history <- "_anwil_score_history\\.csv$"

print(paste("Searching for pattern:", pattern_main))
all_files <- list.files(input_dir, full.names = TRUE)
main_files <- all_files[grepl(pattern_main, all_files)]

print(paste("Found", length(main_files), "main files to plot."))
pdf(output_pdf, width = 10, height = 8)

# Map kernel_type and method_flag to descriptive names
kernel_name <- ifelse(kernel_type == "1", "Gaussian", ifelse(kernel_type == "2", "Tri-cubic", "Unknown"))
method_name <- ifelse(method_flag == "1", "Arithmetic", ifelse(method_flag == "2", "Geometric", "Unknown"))

for (m_file in main_files) {
    dataset_name <- sub("_smoothed.*", "", basename(m_file))
    data_methods <- read_csv(m_file, show_col_types = FALSE)
    
    # Filter metadata files to match the current main file
    base_name <- sub("\\.csv$", "", basename(m_file))
    score_history_file <- all_files[grepl(paste0(base_name, "_anwil_score_history\\.csv$"), all_files)]
    sigma_file <- all_files[grepl(paste0(base_name, pattern_sigma), all_files)]

    print(paste("Processing dataset:", dataset_name))
    print(paste("Main file:", m_file))
    
    # Extract iteration data from score history if available
    stop_iter <- NA
    best_iter <- NA
    best_score <- NA
    if (length(score_history_file) > 0) {
        score_lines <- readLines(score_history_file[1])
        stop_iter <- as.numeric(sub(".*stop_iter=(\\d+).*", "\\1", score_lines[2]))
        best_iter <- as.numeric(sub(".*best_iter=(\\d+).*", "\\1", score_lines[4]))
        best_score <- as.numeric(sub(".*best_score=\\s*([0-9eE+.-]+).*", "\\1", score_lines[5]))
    }

    # 1. Calculate means for centering (x and y) for manifold methods
    x_mean <- mean(data_methods$x_original, na.rm = TRUE)
    y_mean <- mean(data_methods$y_original, na.rm = TRUE)

    # 2. Define strict order for methods
    method_order <- c("Local Sigma (Smooth)", "ManLe", "AmanLe", "ANWIL", "ANWIL (iterative)")

    # 3. Create background dataframes (centered for ManLe/AmanLe, absolute for ANWIL)
    df_background <- bind_rows(
      tibble(x = data_methods$x_original - x_mean, y = data_methods$y_original - y_mean, method = "ManLe"),
      tibble(x = data_methods$x_original - x_mean, y = data_methods$y_original - y_mean, method = "AmanLe"),
      tibble(x = data_methods$x_original,          y = data_methods$y_original,          method = "ANWIL"),
      tibble(x = data_methods$x_original,          y = data_methods$y_original,          method = "ANWIL (iterative)")
    )

    # 4. Create plotting dataframe (Lines for smoothed results)
    df_plot <- bind_rows(
      tibble(x = data_methods$x_manle, y = data_methods$y_manle, method = "ManLe", type = "Main"),
      tibble(x = data_methods$x_amanle, y = data_methods$y_amanle, method = "AmanLe", type = "Main"),
      tibble(x = data_methods$x_anwil, y = data_methods$y_anwil, method = "ANWIL", type = "Main"),
      tibble(x = data_methods$x_anwil_iterative, y = data_methods$y_anwil_iterative, method = "ANWIL (iterative)", type = "Main")
    )

    # Special dataframe for ManLe SVD result
    df_manle_svd <- tibble(
      x = data_methods$x_manle_svd, 
      y = data_methods$y_manle_svd, 
      method = "ManLe", 
      type = "SVD"
    )

    # 5. Handle Sigma visualization (Raw points and Smooth line)
    if (exists("df_sigma_raw_points")) rm(df_sigma_raw_points)
    
    if (length(sigma_file) > 0) {
        data_sigma <- read_csv(sigma_file[1], show_col_types = FALSE)
        df_sigma_long <- tibble(
            x = data_sigma$x, 
            y = data_sigma$local_sigma_smooth, 
            method = "Local Sigma (Smooth)"
        )
        df_sigma_raw_points <- tibble(
            x = data_sigma$x, 
            y = data_sigma$local_sigma_raw,
            method = "Local Sigma (Smooth)"
        )
        df_plot <- bind_rows(df_plot, df_sigma_long)
    }

    # Apply factor levels to force the plot order
    df_plot$method <- factor(df_plot$method, levels = method_order)
    df_background$method <- factor(df_background$method, levels = method_order)
    if (exists("df_sigma_raw_points")) {
        df_sigma_raw_points$method <- factor(df_sigma_raw_points$method, levels = method_order)
    }

    # 6. Generate plot
    p <- ggplot() +
      # Background grey points
      geom_point(data = df_background, aes(x = x, y = y), 
                 color = "gray85", size = 0.7, alpha = 0.8) +
      
      # Main result line for all methods
      geom_line(data = df_plot, aes(x = x, y = y, color = method), 
                linewidth = 0.8) +
      
      # Exclusive SVD line for ManLe (black dashed)
      geom_line(data = df_manle_svd, aes(x = x, y = y), 
                color = "black", linetype = "dashed", linewidth = 0.6) +
      
      # Red points for raw Local Sigma (displayed only in its facet)
      { if(exists("df_sigma_raw_points")) 
          geom_point(data = df_sigma_raw_points, aes(x = x, y = y), 
                     color = "red", size = 0.3, alpha = 0.1) 
      } +
      
      facet_wrap(~ method, scales = "free", ncol = 3) + 
      theme_minimal(base_size = 11) +
      theme(
          legend.position = "none", 
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold")
      ) +
      labs(
        title = paste("Dataset:", dataset_name),
        subtitle = paste(
          "k =", k, "| span =", span, "| iter =", max_iter, "| k_sigma =", k_neighbors_sigma, 
          "| kernel =", kernel_name, "| method =", method_name, 
          "| w_r =", w_r, "| w_e =", w_e, "| w_c =", w_c, "\n",
          "ANWIL iterative: stop_iter =", stop_iter, "| best_iter =", best_iter, "| best_score =", sprintf("%.4f", best_score)
        ),
      )

    print(p)
}

dev.off()
cat("Plots saved at:", output_pdf, "\n")