devtools::load_all()
library(ggplot2)
library(patchwork)

# 1. Simulate Organic Dendritic Neuron
message("Simulating organic branching neuron histology...")
landscape <- simulate_dendritic_neuron(n_grid = 100, noise = 0.15)

# --- VISUALIZATION 1: Raw Organic Channels ---
message("Generating Plot 1: Raw Organic Channels...")

# Neuron/MAP2 - gorgeous emerald/green scale
p_neuron <- ggplot(landscape, aes(x = x, y = y, fill = Neuron)) +
  geom_raster() +
  scale_fill_viridis_c(option = "viridis", name = "MAP2") +
  labs(title = "Neuron (MAP2)", subtitle = "Dendritic arbor with multiple forks") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "gray50"))

# Microglia/Iba1 - hot magma scale showing affinity wrapping the branches
p_iba1 <- ggplot(landscape, aes(x = x, y = y, fill = Microglia)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", name = "Iba1") +
  labs(title = "Microglia (Iba1)", subtitle = "Wrapping active dendritic segments") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "gray50"))

# PNN/WFA - cool mako scale concentrated on soma, excluded from distal forks
p_pnn <- ggplot(landscape, aes(x = x, y = y, fill = PNNs)) +
  geom_raster() +
  scale_fill_viridis_c(option = "mako", name = "WFA") +
  labs(title = "Perineuronal Net (WFA)", subtitle = "Somatic wrapping, distal exclusion") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "gray50"))

# Stitch together and save
combined_channels <- p_neuron + p_iba1 + p_pnn + 
  plot_annotation(
    title = "Organic Branching Neuron Histology Channels",
    subtitle = "Mathematical simulation of branching bifurcations, somatic PNNs, and microglial affinity",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"))
  )

png_channels_path <- "scratch/neuron_channels.png"
ggsave(png_channels_path, combined_channels, width = 14, height = 5.5, dpi = 300)
message("Saved Plot 1 to: ", png_channels_path)


# --- RUN MCMC ESTIMATION ---
message("\nFitting joint spatial models on organic branching structure...")

# 2. Fit Joint Affinity Model for Neuron & Microglia
message("  Fitting Neuron & Microglia (Joint Affinity)...")
fit_aff <- run_joint_mcmc(
  y1 = landscape$Neuron, y2 = landscape$Microglia,
  baseline1 = rep(0, nrow(landscape)), baseline2 = rep(0, nrow(landscape)),
  latent = landscape$LatentDendrite,
  mode = "affinity",
  iterations = 3000,
  init_beta1 = 1, init_beta2 = 1, init_c = 1,
  sd_beta1 = 0.1, sd_beta2 = 0.1, sd_c = 0.1
)
optimal_c_aff <- median(fit_aff$c)
cat("  Optimal Changepoint c (Affinity along branches):", round(optimal_c_aff, 3), "\n")


# --- VISUALIZATION 2: Fitted Dendritic Boundaries ---
message("\nGenerating Plot 2: Fitted Dendritic Boundaries...")

# Add fitted level-set classifications to the data
# Since the boundaries are defined on the LatentDendrite gradient, 
# drawing the contour traces the EXACT bifurcating branch structure!
landscape$Fitted_Boundary <- abs(landscape$LatentDendrite - optimal_c_aff) < 0.12

# Plot 2: Overlaying the fitted organic boundary on top of the neuron map
p_fit <- ggplot(landscape, aes(x = x, y = y)) +
  geom_raster(aes(fill = Neuron)) +
  scale_fill_viridis_c(option = "viridis", name = "MAP2") +
  # Draw estimated organic boundary as a glowing cyan line
  geom_tile(data = subset(landscape, Fitted_Boundary), fill = "cyan", color = "cyan", alpha = 0.9) +
  labs(title = "Fitted Organic Dendritic Boundary", 
       subtitle = "Cyan line: Estimated changepoint boundary dynamically tracking all forks and branches") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray30"))

png_fit_path <- "scratch/fitted_neuron_boundary.png"
ggsave(png_fit_path, p_fit, width = 8, height = 7.5, dpi = 300)
message("Saved Plot 2 to: ", png_fit_path)
