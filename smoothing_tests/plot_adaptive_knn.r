library(ggplot2)

# Read the adaptive KNN data
df <- read.table("adaptive_knn.tsv", header = TRUE, sep = "\t")

# Print the first few rows of the data
print(head(df))

# Plot adaptive KNN smoothing vs LOESS
p <- ggplot(df, aes(x = x)) +
    geom_point(aes(y = y_noisy), alpha = 0.3, color="red") +
    geom_line(aes(y = y_smooth), color="blue", size=0.7) +
    geom_smooth(aes(y = y_noisy), method = "loess", se = FALSE, color = "green", linetype = "dashed") +
    theme_minimal() +
    labs(title="Adaptive KNN smoothing vs LOESS",
         x="x",
         y="y")

# Display the plot
print(p)

# Save the plot as a file
ggsave("adaptive_knn_loess_plot.png", plot = p, width=7, height=5, dpi=400)

# Read curvature and smoothed values data
curvature_df <- read.table("curvature_values_global.tsv", header = TRUE, sep = "\t")
smoothed_df <- read.table("smoothed_values_global.tsv", header = TRUE, sep = "\t")

# Calculate the average curvature for each k
curvature_avg_df <- aggregate(C ~ k, data = curvature_df, FUN = mean)

# Plot the elbow (average curvature vs k)
p_curvature <- ggplot(curvature_avg_df, aes(x = k, y = C)) +
    geom_line(color = "blue", size = 0.7) +
    geom_point(color = "red", size = 1.5) +
    theme_minimal() +
    labs(title = "Average Curvature vs k (Elbow)",
         x = "k",
         y = "Average Curvature (C)")

# Display the curvature plot
print(p_curvature)

# Save the curvature plot as a file
ggsave("curvature_plot.png", plot = p_curvature, width = 7, height = 5, dpi = 400)

# Plot smoothed values for different k with original points
df_noisy <- df[, c("x", "y_noisy")]

p_smoothed <- ggplot() +
    geom_point(data = df_noisy, aes(x = x, y = y_noisy), alpha = 0.3, color = "red") +
    geom_line(data = smoothed_df, aes(x = Point_Index, y = Smoothed_Value, color = as.factor(k)), size = 0.7) +
    theme_minimal() +
    labs(title = "Smoothed values for different k with original points",
         x = "Point Index",
         y = "Smoothed Value",
         color = "k")

# Display the smoothed values plot
print(p_smoothed)

# Save the smoothed values plot as a file
ggsave("smoothed_values_with_originals_plot.png", plot = p_smoothed, width = 7, height = 5, dpi = 400)
