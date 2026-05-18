library(spatialbp)
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(123)

# ==========================================
# 1. Simulate 2D Organic Blobs (Imperfect Rules)
# ==========================================
# Create a 100x100 grid
grid_size <- 100
df <- expand.grid(x = 1:grid_size, y = 1:grid_size)

# Generate an organic "blob" field by summing multiple 2D Gaussians
blob_centers <- data.frame(
  x = c(25, 75, 40, 80),
  y = c(30, 80, 70, 20),
  sigma = c(15, 20, 10, 15),
  weight = c(1.0, 1.2, 0.8, 1.0)
)

df$latent_z <- 0
for (i in 1:nrow(blob_centers)) {
  cx <- blob_centers$x[i]
  cy <- blob_centers$y[i]
  sig <- blob_centers$sigma[i]
  w <- blob_centers$weight[i]
  df$latent_z <- df$latent_z + w * exp(-((df$x - cx)^2 + (df$y - cy)^2) / (2 * sig^2))
}

# Scale latent_z to [0, 100] for easier interpretation
df$latent_z <- (df$latent_z - min(df$latent_z)) / (max(df$latent_z) - min(df$latent_z)) * 100

# True Parameters
true_c <- 40.0       # The boundary is the organic contour where latent_z == 40
true_beta_G <- 2.0   
true_beta_R <- 1.5   
true_beta_B <- 2.5   
noise_sd <- 5.0      # Heavy noise

# Green (Base): Increases INSIDE the blobs (z > c)
hinge_G <- pmax(0, df$latent_z - true_c)
df$Green <- 5.0 + true_beta_G * hinge_G + rnorm(nrow(df), 0, noise_sd)

# Red (Affinity): Increases INSIDE the blobs (z > c)
hinge_R <- pmax(0, df$latent_z - true_c)
df$Red <- 2.0 + true_beta_R * hinge_R + rnorm(nrow(df), 0, noise_sd)

# Blue (Repulsion): Increases OUTSIDE the blobs (z < c)
hinge_B <- pmax(0, -(df$latent_z - true_c))
df$Blue <- 5.0 + true_beta_B * hinge_B + rnorm(nrow(df), 0, noise_sd)

# Prevent negative intensities for plotting
df$Green <- pmax(0, df$Green)
df$Red   <- pmax(0, df$Red)
df$Blue  <- pmax(0, df$Blue)

# ==========================================
# 2. Recover the 2D Organic Boundary & Posterior Probabilities
# ==========================================
message("Running 2D Parameter Recovery and Posterior Probabilities...")

# Helper to compute Log-Likelihood for Independent Models (Null Hypothesis)
fit_indep_G <- run_spatial_mcmc(df$Green, rep(5.0, nrow(df)), df$latent_z, 5000, 1.0, mean(df$latent_z), 0.1, 1.0, 0.0)
fit_indep_R <- run_spatial_mcmc(df$Red, rep(2.0, nrow(df)), df$latent_z, 5000, 1.0, mean(df$latent_z), 0.1, 1.0, 0.0)
fit_indep_B <- run_spatial_mcmc(df$Blue, rep(5.0, nrow(df)), df$latent_z, 5000, 1.0, mean(df$latent_z), 0.1, 1.0, 0.0)

ll_G <- mean(fit_indep_G$log_likelihood)
ll_R <- mean(fit_indep_R$log_likelihood)
ll_B <- mean(fit_indep_B$log_likelihood)

# Helper function for safe posterior probability calculation (Log-Sum-Exp trick)
calc_posterior <- function(ll_aff, ll_rep, ll_null) {
  max_ll <- max(ll_aff, ll_rep, ll_null)
  denom <- exp(ll_aff - max_ll) + exp(ll_rep - max_ll) + exp(ll_null - max_ll)
  p_aff <- exp(ll_aff - max_ll) / denom
  p_rep <- exp(ll_rep - max_ll) / denom
  return(list(aff = p_aff, rep = p_rep))
}

# 1. Test Green-Red (Affinity vs Repulsion)
fit_GR_aff <- run_joint_mcmc(df$Green, df$Red, rep(5.0, nrow(df)), rep(2.0, nrow(df)), df$latent_z, "affinity", 5000, 1, 1, mean(df$latent_z), 0.1, 0.1, 1.0, 0, 0)
fit_GR_rep <- run_joint_mcmc(df$Green, df$Red, rep(5.0, nrow(df)), rep(2.0, nrow(df)), df$latent_z, "repulsion", 5000, 1, 1, mean(df$latent_z), 0.1, 0.1, 1.0, 0, 0)

prob_GR <- calc_posterior(mean(fit_GR_aff$log_likelihood), mean(fit_GR_rep$log_likelihood), ll_G + ll_R)

# 2. Test Green-Blue (Affinity vs Repulsion)
fit_GB_aff <- run_joint_mcmc(df$Green, df$Blue, rep(5.0, nrow(df)), rep(5.0, nrow(df)), df$latent_z, "affinity", 5000, 1, 1, mean(df$latent_z), 0.1, 0.1, 1.0, 0, 0)
fit_GB_rep <- run_joint_mcmc(df$Green, df$Blue, rep(5.0, nrow(df)), rep(5.0, nrow(df)), df$latent_z, "repulsion", 5000, 1, 1, mean(df$latent_z), 0.1, 0.1, 1.0, 0, 0)

prob_GB <- calc_posterior(mean(fit_GB_aff$log_likelihood), mean(fit_GB_rep$log_likelihood), ll_G + ll_B)

# Extract Recovered Boundary from the most likely model (Green-Red Affinity)
burn_in <- 2500
recovered_c <- mean(fit_GR_aff$c[(burn_in+1):5000])

# ==========================================
# 3. 2D Visualization
# ==========================================
# Helper function to plot Raster + Boundaries
plot_raster <- function(data, channel, color_scale, title, subtitle) {
  ggplot(data, aes(x = x, y = y)) +
    geom_raster(aes(fill = !!sym(channel))) +
    scale_fill_gradient(low = "black", high = color_scale) +
    geom_contour(aes(z = latent_z), breaks = true_c, color = "white", linetype = "dashed", linewidth = 1) +
    geom_contour(aes(z = latent_z), breaks = recovered_c, color = "yellow", linewidth = 1.2) +
    coord_equal() +
    labs(title = title, subtitle = subtitle, fill = "Intens") +
    theme_void() + theme(legend.position = "right")
}

# Helper function to plot Topological Contours
plot_topology <- function(data, channel, color_scale) {
  ggplot(data, aes(x = x, y = y)) +
    geom_contour(aes(z = !!sym(channel)), color = color_scale, bins = 10, linewidth = 0.5) +
    geom_contour(aes(z = latent_z), breaks = recovered_c, color = "yellow", linewidth = 1.2) +
    coord_equal() +
    labs(title = paste(channel, "Topology"), subtitle = "Channel Iso-Contours vs Recovered Boundary") +
    theme_void() + theme(panel.background = element_rect(fill="black"))
}

# Row 1: Rasters
p_green <- plot_raster(df, "Green", "green", "Green Channel (Base)", "White=True | Yellow=Recovered")
p_red <- plot_raster(df, "Red", "red", "Red Channel (Affinity)", sprintf("P(Affinity) = %.2f%%", prob_GR$aff * 100))
p_blue <- plot_raster(df, "Blue", "blue", "Blue Channel (Repulsion)", sprintf("P(Repulsion) = %.2f%%", prob_GB$rep * 100))

# Row 2: Topological Contours
t_green <- plot_topology(df, "Green", "green")
t_red <- plot_topology(df, "Red", "red")
t_blue <- plot_topology(df, "Blue", "blue")

# Composite plot
final_plot <- (p_green | p_red | p_blue) / (t_green | t_red | t_blue)
final_plot <- final_plot + plot_annotation(
  title = "2D Parameter Recovery and Posterior Probabilities over Organic Blobs",
  subtitle = "Top Row: Spatial Intensity. Bottom Row: Topological Iso-Contours of the raw data."
)

ggsave("2d_blobs_recovery.png", final_plot, width = 15, height = 10, bg = "black")
message("Saved visualization to 2d_blobs_recovery.png")
