library(spatialbp)
library(ggplot2)

data(meuse, package = "sp")
fit <- fit_spatial_boundary(
  df = meuse,
  y1 = "zinc",
  y2 = "lead",
  latent_col = "dist",
  iterations = 3000
)

p <- plot(fit)
ggsave("C:/Users/bindoffa/.gemini/antigravity/brain/93f67760-4ac4-41dc-91ef-e9758c4812b2/artifacts/meuse_ecotone_boundary.png", p, width = 6, height = 6, dpi = 150)
message("Saved meuse_ecotone_boundary.png successfully!")
