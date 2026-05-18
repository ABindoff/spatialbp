devtools::load_all()
library(ggplot2)

# 1. Simulate cortical tissue
message("Simulating tissue...")
df <- simulate_cortical_tissue(n_grid = 100, noise = 0.25)

# 2. Fit using the ultra-simple high-level wrapper
message("Fitting model using high-level wrapper...")
fit <- fit_spatial_boundary(
  df = df,
  y1 = "Astrocytes",
  y2 = "Microglia",
  latent_type = "radial",
  iterations = 2000
)

# 3. Render and save the S3 plot
message("Generating and saving the S3 overlay plot...")
p <- plot(fit)

png_path <- "scratch/easy_wrapper_boundary.png"
ggsave(png_path, p, width = 8, height = 7.0, dpi = 300)
message("Saved wrapper plot to: ", png_path)
