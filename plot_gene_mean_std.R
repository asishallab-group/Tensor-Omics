# R script to visualize gene expression mean-std relationship
# with LOESS smoothing

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read data
raw_data <- read.table(
  "gene_mean_std_raw.tsv",
  header = TRUE,
  sep = "\t",
  colClasses = c("numeric", "numeric", "character")  # mean, std, gene_id
)

loess_data <- read.table("gene_mean_std_loess.tsv", header=TRUE, sep="\t")

head(raw_data)
head(loess_data)

# Basic statistics
cat("Summary of gene expression means:\n")
print(summary(raw_data$mean))
cat("\nSummary of gene expression standard deviations:\n")
print(summary(raw_data$std))
cat("\nNumber of genes analyzed:", nrow(raw_data), "\n")

# Plot 1: Mean vs STD with LOESS smoothing
p1 <- ggplot(raw_data, aes(x = mean, y = std)) +
  geom_point(alpha = 0.3, size = 0.8, color = "gray50") +
  geom_line(data = loess_data, aes(x = x_query, y = y_span_0.4),
            color = "#E41A1C", linewidth = 1.2, linetype = "solid") +
  geom_line(data = loess_data, aes(x = x_query, y = y_span_0.6),
            color = "#377EB8", linewidth = 1.2, linetype = "dashed") +
  labs(title = "Gene Expression: Mean vs Standard Deviation",
       subtitle = "LOESS smoothing reveals heteroscedasticity pattern",
       x = "Mean Expression Level",
       y = "Standard Deviation",
       caption = paste("n =", nrow(raw_data), "genes")) +
  theme_minimal() +
  theme(legend.position = "none")

# Plot 2: All LOESS spans comparison
loess_long <- loess_data %>%
  pivot_longer(cols = -x_query, names_to = "span", values_to = "std_pred")

p2 <- ggplot() +
  geom_point(data = raw_data, aes(x = mean, y = std),
             alpha = 0.2, size = 0.5, color = "gray70") +
  geom_line(data = loess_long, aes(x = x_query, y = std_pred, color = span),
            linewidth = 1.2) +
  scale_color_manual(
    name = "LOESS span",
    values = c("y_span_0.2" = "#4DAF4A",
               "y_span_0.4" = "#E41A1C",
               "y_span_0.6" = "#377EB8",
               "y_span_0.8" = "#984EA3"),
    labels = c("0.2 (local)", "0.4", "0.6", "0.8 (global)")) +
  labs(title = "LOESS Smoothing with Different Spans",
       subtitle = "Gene expression mean-std relationship",
       x = "Mean Expression Level", y = "Standard Deviation") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save plots
pdf("gene_mean_std_analysis.pdf", width = 10, height = 8)
print(p1)
print(p2)
dev.off()

ggsave("gene_mean_std_main.png", p1, width = 8, height = 6, dpi = 300)
ggsave("gene_mean_std_spans.png", p2, width = 8, height = 6, dpi = 300)

# Calculate correlation
correlation <- cor(raw_data$mean, raw_data$std, use = "complete.obs")
cat("\nCorrelation between mean and std:", round(correlation, 3), "\n")

# Fit linear model for comparison
lm_fit <- lm(std ~ mean, data = raw_data)
cat("\nLinear model summary:\n")
print(summary(lm_fit))

cat("\n=== Analysis Complete ===\n")
cat("Files generated:\n")
cat("• gene_mean_std_raw.tsv - Raw mean/std values\n")
cat("• gene_mean_std_loess.tsv - LOESS smoothed curves\n")
cat("• gene_mean_std_analysis.pdf - PDF with all plots\n")
cat("• gene_mean_std_main.png - Main plot (PNG)\n")
cat("• gene_mean_std_spans.png - Span comparison (PNG)\n")
