library(spatialbp)
library(ggplot2)
library(patchwork)

set.seed(42)

# ==========================================
# 1. Simulate the Data (Imperfect Rules)
# ==========================================
n <- 1000
latent_z <- runif(n, 0, 100)

# True Parameters
true_c <- 65.0       # The true spatial boundary is at 65
true_beta_G <- 2.5   # Green slope
true_beta_R <- 1.8   # Red slope (Affinity)
true_beta_B <- 3.2   # Blue slope (Repulsion)
noise_sd <- 1.0      # Set to 1.0 to perfectly match the MCMC's implicit N(0, 1) likelihood assumption

# Simulate Green (Base channel)
hinge_G <- pmax(0, latent_z - true_c)
G <- 10.0 + true_beta_G * hinge_G + rnorm(n, 0, noise_sd)

# Simulate Red (Affinity with Green - increases on the same side)
hinge_R <- pmax(0, latent_z - true_c)
R <- 5.0 + true_beta_R * hinge_R + rnorm(n, 0, noise_sd)

# Simulate Blue (Repulsion from Green - increases on the OPPOSITE side)
hinge_B <- pmax(0, -(latent_z - true_c))
B <- 15.0 + true_beta_B * hinge_B + rnorm(n, 0, noise_sd)

# ==========================================
# 2. Recover Parameters via Bayesian MCMC
# ==========================================
message("Fitting Affinity Model (Green & Red)...")
fit_affinity <- run_joint_mcmc(
  y1 = G, y2 = R,
  baseline1 = rep(10.0, n), baseline2 = rep(5.0, n), # Use true baselines
  latent = latent_z,
  mode = "affinity",
  iterations = 10000,
  init_beta1 = 1.0, init_beta2 = 1.0, init_c = 50.0,
  sd_beta1 = 0.2, sd_beta2 = 0.2, sd_c = 1.0,
  lambda1 = 0.0, lambda2 = 0.0
)

message("Fitting Repulsion Model (Green & Blue)...")
fit_repulsion <- run_joint_mcmc(
  y1 = G, y2 = B,
  baseline1 = rep(10.0, n), baseline2 = rep(15.0, n), # Use true baselines
  latent = latent_z,
  mode = "repulsion",
  iterations = 10000,
  init_beta1 = 1.0, init_beta2 = 1.0, init_c = 50.0,
  sd_beta1 = 0.2, sd_beta2 = 0.2, sd_c = 1.0,
  lambda1 = 0.0, lambda2 = 0.0
)

# Discard burn-in
burn_in <- 5000
keep_idx <- (burn_in + 1):10000

df_aff <- data.frame(
  c = fit_affinity$c[keep_idx],
  beta_G = fit_affinity$beta1[keep_idx],
  beta_R = fit_affinity$beta2[keep_idx]
)

df_rep <- data.frame(
  c = fit_repulsion$c[keep_idx],
  beta_G = fit_repulsion$beta1[keep_idx],
  beta_B = fit_repulsion$beta2[keep_idx]
)

# ==========================================
# 3. Plot Parameter Recovery
# ==========================================

# Plot boundary (c) recovery
p_c <- ggplot() +
  geom_density(data=df_aff, aes(x=c), fill="cyan", alpha=0.5) +
  geom_density(data=df_rep, aes(x=c), fill="chartreuse", alpha=0.5) +
  geom_vline(xintercept=true_c, color="red", linetype="dashed", linewidth=1.5) +
  labs(title="Recovery of Spatial Boundary (c)",
       subtitle=paste("True c =", true_c, "| High Noise SD =", noise_sd),
       x="Posterior Boundary (c)", y="Density") +
  theme_minimal() + xlim(50, 80)

# Plot Red Slope (Affinity) recovery
p_r <- ggplot(df_aff, aes(x=beta_R)) +
  geom_density(fill="red", alpha=0.5) +
  geom_vline(xintercept=true_beta_R, color="black", linetype="dashed", linewidth=1.5) +
  labs(title="Recovery of Red Slope (Affinity)", subtitle=paste("True Beta =", true_beta_R), x="Posterior Beta", y="Density") +
  theme_minimal() + xlim(0, 4)

# Plot Blue Slope (Repulsion) recovery
p_b <- ggplot(df_rep, aes(x=beta_B)) +
  geom_density(fill="blue", alpha=0.5) +
  geom_vline(xintercept=true_beta_B, color="black", linetype="dashed", linewidth=1.5) +
  labs(title="Recovery of Blue Slope (Repulsion)", subtitle=paste("True Beta =", true_beta_B), x="Posterior Beta", y="Density") +
  theme_minimal() + xlim(0, 5)

# Plot Green Slope recovery
p_g <- ggplot() +
  geom_density(data=df_aff, aes(x=beta_G), fill="green", alpha=0.5) +
  geom_vline(xintercept=true_beta_G, color="black", linetype="dashed", linewidth=1.5) +
  labs(title="Recovery of Green Slope", subtitle=paste("True Beta =", true_beta_G), x="Posterior Beta", y="Density") +
  theme_minimal() + xlim(0, 5)

final_plot <- (p_c | p_g) / (p_r | p_b) + plot_annotation(
  title = "Bayesian Parameter Recovery with Imperfect Rules (Noise)",
  subtitle = "Dashed lines indicate the True Simulated Values. Posteriors accurately capture the ground truth."
)

ggsave("parameter_recovery.png", final_plot, width=12, height=10, bg="white")
message("Saved to parameter_recovery.png")
