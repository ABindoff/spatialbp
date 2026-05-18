library(spatstat.data)
library(spatstat.explore)
library(ggplot2)

data(gorillas)

# Calculate gorilla nest density
nest_dens <- density(gorillas, sigma=20)

# Extract water distance
water_dist <- as.data.frame(gorillas.extra$waterdist)
dens_df <- as.data.frame(nest_dens)

df <- merge(water_dist, dens_df, by=c("x", "y"))
names(df) <- c("x", "y", "Distance_to_Water", "Nest_Density")

p <- ggplot(df, aes(x=Distance_to_Water, y=Nest_Density)) + geom_point(alpha=0.5, color="darkgreen") + theme_minimal()
ggsave("gorillas_water.png", p, width=6, height=5)
