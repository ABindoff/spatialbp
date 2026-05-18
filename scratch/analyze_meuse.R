data(meuse, package = "sp")
cat("Correlation zinc vs dist:", cor(meuse$zinc, meuse$dist), "\n")
cat("Correlation lead vs dist:", cor(meuse$lead, meuse$dist), "\n")
# Print a table of mean zinc by distance bins
dist_bins <- cut(meuse$dist, breaks = seq(0, 1, by = 0.1))
print(tapply(meuse$zinc, dist_bins, mean))
