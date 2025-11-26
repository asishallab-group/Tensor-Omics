# Librerías necesarias
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Lista de valores de k usados en el benchmark
# k_list <- c(50, 80, 110, 150, 200, 300)
k_list <- c(20, 30, 50, 80, 110, 150)

# Leer y combinar todos los archivos
all_data <- lapply(k_list, function(k) {
  fname <- sprintf("knn_smoothing_results_k%d.tsv", k)
  df <- read_tsv(fname, show_col_types = FALSE)
  df$k <- k
  df
}) %>% bind_rows()

# Ordenar por Mean dentro de cada grupo de k
all_data <- all_data %>% 
  group_by(k) %>% 
  arrange(Mean, .by_group = TRUE) %>% 
  ungroup()

# Graficar puntos originales + todas las curvas suavizadas
p <- ggplot() +
  # puntos originales
  geom_point(data = all_data, 
             aes(x = Mean, y = OriginalStd), 
             color = "black", 
             alpha = 0.15, 
             size = 0.5) +
  
  # líneas suavizadas
  geom_line(data = all_data, 
            aes(x = Mean, y = SmoothedStd_KNN_k, 
                color = factor(k), 
                group = k),
            linewidth = 1) +
  
  labs(title = "KNN Smoothing: Smoothed OriginalStd vs Mean",
       x = "Mean",
       y = "Std",
       color = "k") +
  scale_color_discrete(name = "k") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(face="bold"),
    legend.position = "right"
  )

ggsave("knn_smoothing_kallisto.png", p, width = 8, height = 5, dpi = 300)
