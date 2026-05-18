library(spatstat.data)
library(spatstat.explore)
library(spatstat.geom)

data(bei)

# bei is a point pattern of trees.
# bei.extra is a list of pixel images (im objects): elev and grad.

# Let's compute tree density as a pixel image
tree_dens <- density(bei, sigma=10)

# Convert all images to a data frame
elev_df <- as.data.frame(bei.extra$elev)
grad_df <- as.data.frame(bei.extra$grad)
dens_df <- as.data.frame(tree_dens)

# We can merge them by x, y
df <- merge(elev_df, grad_df, by=c("x", "y"), suffixes=c("_elev", "_grad"))
df <- merge(df, dens_df, by=c("x", "y"))
names(df) <- c("x", "y", "elev", "grad", "tree_density")

print(head(df))
print(summary(df))

# Is there a correlation?
print(cor(df$elev, df$tree_density, use="complete.obs"))
print(cor(df$grad, df$tree_density, use="complete.obs"))

# Let's save a quick plot to see it!
library(ggplot2)
p <- ggplot(df, aes(x=x, y=y, fill=tree_density)) + geom_raster() + scale_fill_viridis_c() + theme_minimal()
ggsave("bei_test.png", p, width=6, height=4)
