library(ggplot2)
library(spatstat.data)
library(spatstat.explore)
devtools::load_all(".")

data(sporophores)
fungi_dens <- density(sporophores, sigma = 10, eps = 1)
df <- as.data.frame(fungi_dens)
names(df) <- c("x", "y", "Fungus_Density")
df$Distance_To_Tree <- sqrt(df$x^2 + df$y^2)
df <- df[df$Distance_To_Tree <= 50, ]

fit_fungi <- fit_spatial_boundary(
  df = df,
  y1 = "Fungus_Density",
  latent_col = "Distance_To_Tree",
  iterations = 3000
)

cat("Estimated fungal boundary radius:", fit_fungi$optimal_c, "meters from the tree.\n")

p_diag <- plot_diagnostics(fit_fungi, burn_in = 0.2)
ggsave("fungi_diagnostics.png", p_diag[[1]], width=6, height=4) # save trace plot
