library(readr)
library(ggplot2)

df <- read_tsv("synthetic_points.tsv")

df$loess_005 <- loess(y_noisy ~ x, data=df, span=0.05)$fitted
df$loess_02  <- loess(y_noisy ~ x, data=df, span=0.20)$fitted
df$loess_04  <- loess(y_noisy ~ x, data=df, span=0.40)$fitted
df$loess_06  <- loess(y_noisy ~ x, data=df, span=0.60)$fitted

ggplot(df, aes(x=x)) +
  geom_point(aes(y=y_noisy), size=0.5, alpha=0.4) +
  geom_line(aes(y=loess_005, color="span=0.05"), linewidth=1) +
  geom_line(aes(y=loess_02,  color="span=0.20"), linewidth=1) +
  geom_line(aes(y=loess_04,  color="span=0.40"), linewidth=1) +
  geom_line(aes(y=loess_06,  color="span=0.60"), linewidth=1) +
  scale_color_manual(values=c(
    "span=0.05"="red",
    "span=0.20"="orange",
    "span=0.40"="green",
    "span=0.60"="blue"
  )) +
  labs(title="LOESS aplicado a los datos EXACTOS de Fortran",
       y="y_noisy",
       color="Smoothing") +
  theme_minimal()
