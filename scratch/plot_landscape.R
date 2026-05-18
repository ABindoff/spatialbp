source("../R/simulate_landscape.R")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr", repos = "https://cloud.r-project.org")

library(ggplot2)
library(tidyr)

# Generate the data
set.seed(42)
landscape <- simulate_landscape(n_grid = 150, noise = 0.3)

# Reshape for ggplot
landscape_long <- pivot_longer(landscape, 
                               cols = c("Red", "Green", "Blue"),
                               names_to = "Color", 
                               values_to = "Intensity")

# Create the plot
p <- ggplot(landscape_long, aes(x = x, y = y, fill = Intensity)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") +
  facet_wrap(~ Color) +
  theme_minimal() +
  coord_fixed() +
  labs(title = "Simulated Alien Landscape (Satellite Sensor Data)",
       subtitle = "Notice Green and Blue are strictly disjoint, while Red encompasses both.",
       x = "Longitude", y = "Latitude")

ggsave("../landscape_preview.png", p, width = 10, height = 4, bg = "white")
cat("Plot saved to landscape_preview.png\n")
