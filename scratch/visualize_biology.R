devtools::load_all()
library(ggplot2)
library(patchwork)

# 1. Simulate Cortical Histological Tissue
message("Simulating cortical histological tissue (Astrocytes, Microglia, PNNs)...")
landscape <- simulate_cortical_tissue(n_grid = 100, noise = 0.2)

# --- VISUALIZATION 1: Raw Biological Channels ---
message("Generating Plot 1: Raw Histological Channels...")

# Astrocyte channel (GFAP) - warm copper/orange scale
p_gfap <- ggplot(landscape, aes(x = x, y = y, fill = Astrocytes)) +
  geom_raster() +
  scale_fill_viridis_c(option = "inferno", name = "GFAP") +
  labs(title = "Astrocytes (GFAP)", subtitle = "Reactive central gliosis") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "gray50"))

# Microglia channel (Iba1) - hot magma scale
p_iba1 <- ggplot(landscape, aes(x = x, y = y, fill = Microglia)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", name = "Iba1") +
  labs(title = "Microglia (Iba1)", subtitle = "Recruited to gliosis hotspot") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "gray50"))

# PNN channel (WFA) - cool viridis scale
p_pnn <- ggplot(landscape, aes(x = x, y = y, fill = PNNs)) +
  geom_raster() +
  scale_fill_viridis_c(option = "mako", name = "WFA") +
  labs(title = "Perineuronal Networks (WFA)", subtitle = "Excluded & degraded in center") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "gray50"))

# Stitch together and save
combined_channels <- p_gfap + p_iba1 + p_pnn + 
  plot_annotation(
    title = "Simulated 2D Cortical Histology Channels",
    subtitle = "Simulating neuroinflammatory gliosis and perineuronal network degradation",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"))
  )

png_channels_path <- "scratch/cortical_tissue_channels.png"
ggsave(png_channels_path, combined_channels, width = 14, height = 5.5, dpi = 300)
message("Saved Plot 1 to: ", png_channels_path)


# --- RUN MCMC ESTIMATION ---
message("\nFitting joint spatial models via Rust MCMC...")

# 2. Fit Joint Affinity Model for Astrocytes & Microglia
# They are both active in the center (Z > c)
message("  Fitting Astrocytes & Microglia (Joint Affinity)...")
fit_aff <- run_joint_mcmc(
  y1 = landscape$Astrocytes, y2 = landscape$Microglia,
  baseline1 = rep(0, nrow(landscape)), baseline2 = rep(0, nrow(landscape)),
  latent = landscape$LatentInflame,
  mode = "affinity",
  iterations = 3000,
  init_beta1 = 1, init_beta2 = 1, init_c = 0,
  sd_beta1 = 0.1, sd_beta2 = 0.1, sd_c = 0.1
)
optimal_c_aff <- median(fit_aff$c)
cat("  Optimal Changepoint c (Affinity):", round(optimal_c_aff, 3), "\n")

# 3. Fit Joint Repulsion Model for Microglia & PNNs
# Microglia is active in center (Z > c), PNNs active on periphery (Z < c)
message("  Fitting Microglia & PNNs (Joint Repulsion)...")
fit_rep <- run_joint_mcmc(
  y1 = landscape$Microglia, y2 = landscape$PNNs,
  baseline1 = rep(0, nrow(landscape)), baseline2 = rep(0, nrow(landscape)),
  latent = landscape$LatentInflame,
  mode = "repulsion",
  iterations = 3000,
  init_beta1 = 1, init_beta2 = 1, init_c = 0,
  sd_beta1 = 0.1, sd_beta2 = 0.1, sd_c = 0.1
)
optimal_c_rep <- median(fit_rep$c)
cat("  Optimal Changepoint c (Repulsion):", round(optimal_c_rep, 3), "\n")


# --- VISUALIZATION 2: Fitted Spatial Boundaries ---
message("\nGenerating Plot 2: Fitted Model Boundaries...")

# Add fitted level-set classifications to the data
landscape$Fitted_Aff_Boundary <- abs(landscape$LatentInflame - optimal_c_aff) < 0.05
landscape$Fitted_Rep_Boundary <- abs(landscape$LatentInflame - optimal_c_rep) < 0.05

# Plot 2a: Astrocytes + Microglia with fitted joint boundary
p_fit_aff <- ggplot(landscape, aes(x = x, y = y)) +
  geom_raster(aes(fill = Microglia)) +
  scale_fill_viridis_c(option = "magma", name = "Iba1") +
  geom_tile(data = subset(landscape, Fitted_Aff_Boundary), fill = "cyan", color = "cyan", alpha = 0.8) +
  labs(title = "Joint Affinity Boundary", 
       subtitle = "Cyan line: Shared boundary between Astrocytes & Microglia") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "gray30"))

# Plot 2b: Microglia + PNNs with fitted joint boundary
p_fit_rep <- ggplot(landscape, aes(x = x, y = y)) +
  geom_raster(aes(fill = PNNs)) +
  scale_fill_viridis_c(option = "mako", name = "WFA") +
  geom_tile(data = subset(landscape, Fitted_Rep_Boundary), fill = "chartreuse", color = "chartreuse", alpha = 0.8) +
  labs(title = "Joint Repulsion Boundary", 
       subtitle = "Green line: Boundary segregating Microglia and PNNs") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "gray30"))

# Stitch together and save
combined_fits <- p_fit_aff + p_fit_rep +
  plot_annotation(
    title = "Fitted Joint Spatial boundaries (Rust MCMC)",
    subtitle = "Tying boundary contours across channels dynamically isolates affinity and repulsion zones",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"))
  )

png_fits_path <- "scratch/fitted_spatial_boundaries.png"
ggsave(png_fits_path, combined_fits, width = 12, height = 6.5, dpi = 300)
message("Saved Plot 2 to: ", png_fits_path)
