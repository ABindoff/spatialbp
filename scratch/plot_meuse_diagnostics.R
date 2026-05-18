library(spatialbp)
library(ggplot2)

data(meuse, package = "sp")

fit <- fit_spatial_boundary(
  df = meuse,
  y1 = "zinc",
  y2 = "lead",
  latent_col = "dist",
  iterations = 10000
)

# Render diagnostics
png("C:/Users/bindoffa/.gemini/antigravity/brain/93f67760-4ac4-41dc-91ef-e9758c4812b2/artifacts/meuse_diagnostics.png", width = 800, height = 800)
plot_diagnostics(fit)
dev.off()

# Render scatter plot Zinc vs Dist
p_grad <- ggplot(meuse, aes(x = dist, y = zinc, color = lead)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis_c(option = "magma", name = "Lead (ppm)") +
  geom_vline(xintercept = fit$optimal_c, color = "cyan", linetype = "dashed", size = 1.2) +
  labs(
    title = "Heavy Metal Concentration vs. Distance to River Bank",
    subtitle = paste("Estimated Ecotone Boundary at dist =", round(fit$optimal_c, 3)),
    x = "Normalized Distance to River Bank",
    y = "Zinc Concentration (ppm)"
  ) +
  theme_minimal()

ggsave("C:/Users/bindoffa/.gemini/antigravity/brain/93f67760-4ac4-41dc-91ef-e9758c4812b2/artifacts/meuse_scatter.png", p_grad, width = 6, height = 5)
