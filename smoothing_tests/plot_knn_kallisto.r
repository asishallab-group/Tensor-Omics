# Load required libraries
library(ggplot2)

# Read the file
data <- read.table("adaptive_knn_smoothing_results.tsv", header = TRUE, sep = "\t")

# Create the plot
plot <- ggplot(data, aes(x = Mean)) +
  # Original points
  geom_point(aes(y = OriginalStd), color = "blue", alpha = 0.5, size = 1.5) +
  # Smoothed curve
  geom_line(aes(y = SmoothedStd_AdaptiveKNN), color = "red", size = 1) +
  labs(
    title = "Adaptive KNN Smoothing Results",
    x = "Gene Mean",
    y = "Standard Deviation",
    caption = "Blue: Original Std, Red: Smoothed Std"
  ) +
  theme_minimal()

# Display the plot
print(plot)

# Save the plot as a file
ggsave("adaptive_knn_plot.png", plot, width = 8, height = 6)