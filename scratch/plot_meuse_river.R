library(spatialbp)
library(ggplot2)

data(meuse, package = "sp")
data(meuse.riv, package = "sp")

fit <- fit_spatial_boundary(
  df = meuse,
  y1 = "zinc",
  y2 = "lead",
  latent_col = "dist",
  iterations = 3000
)

# Convert river coordinates to a data frame
riv_df <- as.data.frame(meuse.riv)
colnames(riv_df) <- c("x", "y")

# Get default S3 plot
p <- plot(fit)

# Append river polygon layer
p_with_river <- p +
  geom_polygon(data = riv_df, aes(x = x, y = y), fill = "aliceblue", color = "dodgerblue", alpha = 0.6) +
  # Restrict coord limits to the range of the meuse coordinates so the river doesn't zoom out the plot!
  coord_equal(xlim = range(meuse$x), ylim = range(meuse$y))

ggsave("C:/Users/bindoffa/.gemini/antigravity/brain/93f67760-4ac4-41dc-91ef-e9758c4812b2/artifacts/meuse_river_boundary.png", p_with_river, width = 7, height = 7, dpi = 150)
message("Saved meuse_river_boundary.png successfully!")
