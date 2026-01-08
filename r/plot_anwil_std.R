# Load necessary libraries
library(ggplot2)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_anwil_std.R <k_value> <max_iter>")
}

# Extract k and max_iter from arguments
k <- args[1]
max_iter <- args[2]

# Define the input and output directories
input_dir <- "results/data"
output_dir <- "results/plots"
output_pdf <- file.path(output_dir, paste0("anwil_std_plots_k", k, "_iter", max_iter, ".pdf"))

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get the list of CSV files matching "anwil_std"
csv_files <- list.files(input_dir, pattern = paste0("k", k, "_iter", max_iter, "_anwil_std.*\\.csv$"), full.names = TRUE)

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
  p <- ggplot(data, aes(x = x, y = local_sigma)) +
    geom_line(color = "blue") +
    ggtitle(paste("Dataset:", dataset_name, "| k=", k, ", iter=", max_iter)) +
    xlab("X") +
    ylab("Local Sigma") +
    theme_minimal()

  # Print the plot to the PDF device
  print(p)
}

# Close the PDF device
dev.off()

cat("Plots saved to:", output_pdf, "\n")
