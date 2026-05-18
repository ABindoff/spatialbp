library(spatialbp)
library(ggplot2)
library(spatstat.data)
library(spatstat.explore)

data(bei)

# Compute tree density as a pixel image
tree_dens <- density(bei, sigma=15)

# Convert all images to a data frame
elev_df <- as.data.frame(bei.extra$elev)
dens_df <- as.data.frame(tree_dens)

df <- merge(elev_df, dens_df, by=c("x", "y"))
names(df) <- c("x", "y", "Elevation", "Tree_Density")

p <- ggplot(df, aes(x=Elevation, y=Tree_Density)) + geom_point(alpha=0.1) + theme_minimal()
ggsave("bei_elev_scatter.png", p, width=6, height=5)
