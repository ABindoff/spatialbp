library(spatialbp)
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(42)

# ==========================================
# 1. Simulate Competitive Co-Dependence
# ==========================================
# Green is a large structural blob in the center.
# Red and Blue both require Green (Radial Affinity).
# But Red and Blue compete, splitting the territory left/right (Horizontal Repulsion).

grid_size <- 100
df <- expand.grid(x = 1:grid_size, y = 1:grid_size)

center_x <- 50
center_y <- 50
df$radial_z <- sqrt((df$x - center_x)^2 + (df$y - center_y)^2)
df$horizontal_z <- df$x

noise_sd <- 2.0

# Green: High in the center (radial distance < 35)
hinge_G <- pmax(0, -(df$radial_z - 35))
df$Green <- 5.0 + 3.0 * hinge_G + rnorm(nrow(df), 0, noise_sd)

# Red: High in the center, BUT restricted to the left (x < 50)
hinge_R_rad <- pmax(0, -(df$radial_z - 35))
hinge_R_hor <- pmax(0, -(df$horizontal_z - 50))
df$Red <- 2.0 + 2.0 * hinge_R_rad + 2.0 * hinge_R_hor + rnorm(nrow(df), 0, noise_sd)

# Blue: High in the center, BUT restricted to the right (x > 50)
hinge_B_rad <- pmax(0, -(df$radial_z - 35))
hinge_B_hor <- pmax(0, (df$horizontal_z - 50))
df$Blue <- 2.0 + 2.0 * hinge_B_rad + 2.0 * hinge_B_hor + rnorm(nrow(df), 0, noise_sd)

# Prevent negative intensities
df$Green <- pmax(0, df$Green)
df$Red   <- pmax(0, df$Red)
df$Blue  <- pmax(0, df$Blue)

# Scale data for discovery
df$Green_scaled <- scale(df$Green)
df$Red_scaled <- scale(df$Red)
df$Blue_scaled <- scale(df$Blue)

stains <- list(Green = df$Green_scaled, Red = df$Red_scaled, Blue = df$Blue_scaled)

# ==========================================
# 2. Discover Relationships
# ==========================================
message("Testing Radial Axis (Core vs Periphery)...")
res_radial <- discover_spatial_relationships(stains, latent_field = df$radial_z, iterations = 2000)

message("Testing Horizontal Axis (Left vs Right)...")
res_horizontal <- discover_spatial_relationships(stains, latent_field = df$horizontal_z, iterations = 2000)

# ==========================================
# 3. Build Interaction Graph Plot
# ==========================================
# We manually create a node/edge graph visualization
nodes <- data.frame(
  name = c("Green", "Red", "Blue"),
  x = c(50, 30, 70),
  y = c(80, 20, 20),
  color = c("#00FF00", "#FF0000", "#0000FF")
)

edges <- data.frame(
  from = c("Green", "Green", "Red"),
  to = c("Red", "Blue", "Blue"),
  type = c("Affinity", "Affinity", "Repulsion"),
  x = c(50, 50, 30),
  y = c(80, 80, 20),
  xend = c(30, 70, 70),
  yend = c(20, 20, 20)
)

p_graph <- ggplot() +
  geom_segment(data = edges %>% filter(type == "Affinity"), 
               aes(x=x, y=y, xend=xend, yend=yend), color="white", linewidth=2) +
  geom_segment(data = edges %>% filter(type == "Repulsion"), 
               aes(x=x, y=y, xend=xend, yend=yend), color="red", linewidth=2, linetype="dashed") +
  geom_point(data = nodes, aes(x=x, y=y, color=color), size=20) +
  geom_text(data = nodes, aes(x=x, y=y, label=name), color="black", fontface="bold") +
  annotate("text", x = 40, y = 50, label = "Radial Affinity", color="white", angle=60) +
  annotate("text", x = 60, y = 50, label = "Radial Affinity", color="white", angle=-60) +
  annotate("text", x = 50, y = 15, label = "Horizontal Repulsion", color="red") +
  scale_color_identity() +
  theme_void() +
  theme(panel.background = element_rect(fill="black")) +
  labs(title = "Interaction Web: Competitive Co-Dependence")

# ==========================================
# 4. 2D Visualization
# ==========================================
plot_raster <- function(data, channel, color_scale, title) {
  ggplot(data, aes(x = x, y = y)) +
    geom_raster(aes(fill = !!sym(channel))) +
    scale_fill_gradient(low = "black", high = color_scale) +
    coord_equal() +
    labs(title = title, fill = "Intens") +
    theme_void() + theme(legend.position = "none")
}

p_green <- plot_raster(df, "Green", "green", "Green (Structural Base)")
p_red <- plot_raster(df, "Red", "red", "Red (Competitor A)")
p_blue <- plot_raster(df, "Blue", "blue", "Blue (Competitor B)")

final_plot <- (p_green | p_red | p_blue) / p_graph
final_plot <- final_plot + plot_annotation(
  title = "Competitive Co-Dependence in Spatial Communities",
  subtitle = "Red & Blue both require Green's radial territory (Affinity), but repel each other horizontally."
)

ggsave("competitive_codependence.png", final_plot, width = 12, height = 8, bg = "black")
message("Saved visualization to competitive_codependence.png")
