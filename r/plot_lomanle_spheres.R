#!/usr/bin/env Rscript
library(ggplot2)
library(ggforce)
library(viridis)

# --- Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: Rscript plot_report.R <file.csv> <k_min> <g_threshold> <o_max> <o_min>")

# --- Data Loading ---
df <- read.csv(args[1])
# Prepare output PDF path
pdf_out <- file.path("results/plots", paste0(gsub(".csv", "", basename(args[1])), "_5page_report.pdf"))
pdf(pdf_out, width = 11, height = 8.5)

# --- PAGE 1: ALL ADAPTIVE SPHERES (Search Space) ---
# Visualizes the local neighborhood size (radius) found by the adaptive k-search
p1 <- ggplot(df) +
  geom_circle(aes(x0 = x, y0 = y, r = radius), color = "blue", alpha = 0.05, linewidth = 0.05) +
  geom_point(aes(x = x, y = y), size = 0.4, color = "black") +
  coord_fixed() + 
  theme_minimal() +
  labs(title = "All Adaptive Spheres", 
       subtitle = "Visualizing the local neighborhood size (radius) for every data point.")
print(p1)

# --- PAGE 2: DENSITY DISTRIBUTION ---
# Warm colors represent higher data concentration (rho)
p2 <- ggplot(df) +
  geom_point(aes(x = x, y = y, color = density), size = 1.2) +
  scale_color_viridis(option = "magma") +
  coord_fixed() + 
  theme_minimal() +
  labs(title = "Point Density Distribution", 
       subtitle = "Warmer colors indicate regions with higher data concentration.")
print(p2)

# --- PAGE 3: GLOBAL TANGENT FIELD (All PC1s) ---
# This visualizes the SVD results for EVERY point, centered at p_i
p3 <- ggplot(df) +
  # Background points (faint)
  geom_point(aes(x = x, y = y), color = "grey90", size = 0.2, alpha = 0.5) +
  
  # PC1 vectors for all points, centered at (x,y)
  geom_segment(aes(x = x - v1_x * s1, 
                   y = y - v1_y * s1, 
                   xend = x + v1_x * s1, 
                   yend = y + v1_y * s1,
                   color = density), 
               linewidth = 0.2, alpha = 0.4) +
  
  scale_color_viridis(option = "magma", name = "Density") +
  coord_fixed() + 
  theme_minimal() +
  labs(title = "Global Tangent Field (All PC1s)",
       subtitle = "Local linear approximation (SVD) centered at each point p_i.")
print(p3)

# --- PAGE 4: ATLAS SELECTION (Anchors with Spheres) ---
# Highlights the points selected as anchors for the Atlas
p4 <- ggplot(df) +
  geom_point(aes(x = x, y = y), color = "grey90", size = 0.3) +
  geom_circle(data = subset(df, anchor == 1),
              aes(x0 = x, y0 = y, r = radius, fill = density), 
              color = "darkred", alpha = 0.2, linewidth = 0.3) +
  geom_point(data = subset(df, anchor == 1), aes(x = x, y = y), color = "black", size = 1) +
  scale_fill_viridis(option = "magma") +
  coord_fixed() + 
  theme_minimal() +
  labs(title = "Atlas Selection (Anchors)")
print(p4)

# --- PAGE 5: ANCHOR SKELETON (PC1 Structure) ---
# Final skeletal structure of the manifold using only Anchor points
p5 <- ggplot(df) +
  geom_point(aes(x = x, y = y), color = "grey90", size = 0.3) +
  # Highlighted segments for anchors only
  geom_segment(data = subset(df, anchor == 1),
               aes(x = x - v1_x * s1, y = y - v1_y * s1, 
                   xend = x + v1_x * s1, yend = y + v1_y * s1),
               color = "red", linewidth = 1) +
  geom_point(data = subset(df, anchor == 1), aes(x = x, y = y), color = "black", size = 1.5) +
  coord_fixed() + 
  theme_minimal() +
  labs(title = "Final Anchor Skeleton")
print(p5)

dev.off()