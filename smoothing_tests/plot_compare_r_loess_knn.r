library(ggplot2)
library(dplyr)
library(readr)

k_list <- c(25,50,100,200)

all <- lapply(k_list, function(k){
  df <- read_tsv(sprintf("knn_sim_k%d.tsv", k), show_col_types = FALSE)
  df$k <- k
  df
}) |> bind_rows()

p <- ggplot(all, aes(x = Mean_x)) +
  geom_point(aes(y = Y_noisy), alpha = 0.2, size = 0.5) +
  geom_line(aes(y = Y_smoothed, color = factor(k)), linewidth = 1) +
  theme_minimal()

ggsave("knn_smoothing_k.png", p, width = 8, height = 5, dpi = 300)
