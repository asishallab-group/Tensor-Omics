# Load necessary libraries
library(ggplot2)

print("Starting to generate plots for anwil_std...")
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript plot_anwil_std.R <k_neighbors> <span> <max_iter> <kernel_type> <k_neighbors_sigma> <method_flag> <w_r> <w_e> <w_c>")
}

# Extract k_neighbors, span, and n_iters_max from arguments
k <- args[1]
span <- args[2]
max_iter <- args[3]
kernel_type <- args[4]
k_neighbors_sigma <- args[5]
method_flag <- args[6]
w_r <- sprintf("%.2f", as.numeric(args[7]))
w_e <- sprintf("%.2f", as.numeric(args[8]))
w_c <- sprintf("%.2f", as.numeric(args[9]))

# Define the input and output directories
input_dir <- "results/data/2d"
output_dir <- "results/plots"
output_pdf <- file.path(output_dir, paste0("anwil_std_plots_k", k, "_iter", max_iter, "_span", format(span, nsmall = 2), "_ksigma", k_neighbors_sigma, "_kernel", kernel_type, "_method", method_flag, "_wr", w_r, "_we", w_e, "_wc", w_c, ".pdf"))

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get the list of CSV files matching "anwil_std"
csv_files <- list.files(input_dir, pattern = paste0("k",k,"_iter",max_iter,"_span",format(span, nsmall = 2),"_ksigma", k_neighbors_sigma, "_kernel", kernel_type, "_method", method_flag, "_wr", w_r, "_we", w_e, "_wc", w_c, "_anwil_std.*\\.csv$"), full.names = TRUE)

# Open a PDF device to save all plots together
pdf(output_pdf)

# Loop through each CSV file and generate a plot
for (csv_file in csv_files) {
  # Read the data
  data <- read.csv(csv_file)

  # Extract the base name of the file (without extension)
  base_name <- tools::file_path_sans_ext(basename(csv_file))

  # Extract the dataset name from the filename
  dataset_name <- sub("_smoothed.*", "", basename(csv_file))

  # Generate the plot
  p <- ggplot(data) +
    geom_point(aes(x = x, y = local_sigma_raw), color = "red", size = 1.5) +
    geom_line(aes(x = x, y = local_sigma_smooth), color = "blue") +
    ggtitle(paste("Dataset:", dataset_name, "| k=", k, ", iter=", max_iter, ", span=", span, ", k_sigma=", k_neighbors_sigma, ", kernel=", kernel_type)) +
    xlab("X") +
    ylab("Local Sigma") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10))

  # Print the plot to the PDF device
  print(p)
}

# Close the PDF device
dev.off()

cat("Plots saved to:", output_pdf, "\n")
