library(spatialbp)
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(42)

# ==========================================
# 1. Simulate 2D Spatial Grid (Imperfect Rules)
# ==========================================
# Create a 100x100 grid
grid_size <- 100
df <- expand.grid(x = 1:grid_size, y = 1:grid_size)

# Define the Latent Variable (Z) as the distance from the center
center_x <- 50
center_y <- 50
df$latent_z <- sqrt((df$x - center_x)^2 + (df$y - center_y)^2)

# True Parameters
true_c <- 25.0       # The true boundary is a circle of radius 25
true_beta_G <- 2.0   
true_beta_R <- 1.5   
true_beta_B <- 2.5   
noise_sd <- 2.0      # Imperfect rules (noise)

# Green (Base): Increases outside the boundary
hinge_G <- pmax(0, df$latent_z - true_c)
df$Green <- 5.0 + true_beta_G * hinge_G + rnorm(nrow(df), 0, noise_sd)

# Red (Affinity): Increases outside the boundary
hinge_R <- pmax(0, df$latent_z - true_c)
df$Red <- 2.0 + true_beta_R * hinge_R + rnorm(nrow(df), 0, noise_sd)

# Blue (Repulsion): Increases INSIDE the boundary
hinge_B <- pmax(0, -(df$latent_z - true_c))
df$Blue <- 5.0 + true_beta_B * hinge_B + rnorm(nrow(df), 0, noise_sd)

# Prevent negative intensities for plotting
df$Green <- pmax(0, df$Green)
df$Red   <- pmax(0, df$Red)
df$Blue  <- pmax(0, df$Blue)

# ==========================================
# 2. Recover the 2D Boundary via MCMC
# ==========================================
message("Running 2D Parameter Recovery...")

# We use the Rust MCMC backend to recover the boundary using Green & Red (Affinity)
fit_affinity <- run_joint_mcmc(
  y1 = df$Green, y2 = df$Red,
  baseline1 = rep(5.0, nrow(df)), baseline2 = rep(2.0, nrow(df)),
  latent = df$latent_z,
  mode = "affinity",
  iterations = 5000,
  init_beta1 = 1.0, init_beta2 = 1.0, init_c = mean(df$latent_z),
  sd_beta1 = 0.1, sd_beta2 = 0.1, sd_c = 1.0,
  lambda1 = 0.0, lambda2 = 0.0
)

# Extract Posterior Mean of c (discarding burn-in)
burn_in <- 2500
recovered_c <- mean(fit_affinity$c[(burn_in+1):5000])
ci_lower_c <- quantile(fit_affinity$c[(burn_in+1):5000], 0.025)
ci_upper_c <- quantile(fit_affinity$c[(burn_in+1):5000], 0.975)

message(sprintf("True Boundary: %.2f | Recovered Boundary: %.2f (95%% CI: %.2f - %.2f)", 
                true_c, recovered_c, ci_lower_c, ci_upper_c))

# ==========================================
# 3. 2D Visualization
# ==========================================
# Helper function to plot a channel
plot_channel <- function(data, channel, color_scale, title) {
  ggplot(data, aes(x = x, y = y, fill = !!sym(channel))) +
    geom_raster() +
    scale_fill_gradient(low = "black", high = color_scale) +
    # Draw the True Boundary (Dashed White)
    annotate("path",
             x = center_x + true_c * cos(seq(0, 2*pi, length.out=100)),
             y = center_y + true_c * sin(seq(0, 2*pi, length.out=100)),
             color = "white", linetype = "dashed", linewidth = 1) +
    # Draw the Recovered MCMC Boundary (Solid Yellow)
    annotate("path",
             x = center_x + recovered_c * cos(seq(0, 2*pi, length.out=100)),
             y = center_y + recovered_c * sin(seq(0, 2*pi, length.out=100)),
             color = "yellow", linewidth = 1.2) +
    coord_equal() +
    labs(title = title, subtitle = "White=True Boundary | Yellow=Recovered", fill = "Intensity") +
    theme_void() +
    theme(legend.position = "right")
}

p_green <- plot_channel(df, "Green", "green", "Green Channel (Base)")
p_red <- plot_channel(df, "Red", "red", "Red Channel (Affinity with Green)")
p_blue <- plot_channel(df, "Blue", "blue", "Blue Channel (Repulsion from Green)")

# Composite plot
final_plot <- p_green | p_red | p_blue
final_plot <- final_plot + plot_annotation(
  title = "2D Spatial Parameter Recovery of Multi-Species Communities",
  subtitle = "Demonstrating how imperfect rules (noise) map to 2D space and how the MCMC engine recovers the exact geographical boundary."
)

ggsave("2d_parameter_recovery.png", final_plot, width = 15, height = 5, bg = "black")
message("Saved visualization to 2d_parameter_recovery.png")
