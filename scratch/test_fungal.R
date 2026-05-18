library(spatialbp)
library(ggplot2)

simulate_fungal_symbiosis <- function(n_grid = 100, noise = 0.15) {
  # Create coordinate grid
  df <- expand.grid(x = seq(0, 1, length.out = n_grid), 
                    y = seq(0, 1, length.out = n_grid))
  
  # Latent gradient: Distance from an ancient Birch tree cluster at (0.3, 0.7)
  # Actually, let's use a multi-modal grove
  dist1 <- sqrt((df$x - 0.3)^2 + (df$y - 0.7)^2)
  dist2 <- sqrt((df$x - 0.7)^2 + (df$y - 0.4)^2)
  # Z is proximity (inverted distance)
  df$Root_Proximity <- 1.0 - pmin(dist1, dist2)
  
  # True ecological boundary: Symbiosis breaks down if proximity < 0.65 (distance > 0.35)
  c_true <- 0.65
  
  # Birch Tree Density (smoothly decays with distance)
  df$Birch_Density <- pmax(0, df$Root_Proximity - 0.4) * 100 + rnorm(nrow(df), 0, noise * 20)
  df$Birch_Density <- pmax(0, df$Birch_Density)
  
  # Mycorrhizal Fungi Spores (Sharp spatial boundary constraint!)
  # They require a minimum root proximity to survive
  fungus_base <- 10
  fungus_active <- ifelse(df$Root_Proximity > c_true, 
                          (df$Root_Proximity - c_true) * 300, 
                          0)
  df$Fungus_Spores <- fungus_base + fungus_active + rnorm(nrow(df), 0, noise * 10)
  df$Fungus_Spores <- pmax(0, df$Fungus_Spores)
  
  return(df)
}

# Generate data
set.seed(42)
df_forest <- simulate_fungal_symbiosis(n_grid = 80)

# Fit boundary
fit <- fit_spatial_boundary(
  df = df_forest,
  y1 = "Birch_Density",
  y2 = "Fungus_Spores",
  latent_col = "Root_Proximity",
  iterations = 3000
)

cat("True c:", 0.65, "\n")
cat("Estimated c:", fit$optimal_c, "\n")

# Plot boundary
p <- plot(fit)
ggsave("fungal_symbiosis_map.png", p, width=6, height=5)

p_diag <- plot_diagnostics(fit)
ggsave("fungal_diagnostics.png", p_diag[[1]], width=6, height=4) # just save the trace
