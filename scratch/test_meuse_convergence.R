library(spatialbp)

data(meuse, package = "sp")

# Test with 10,000 iterations
fit <- fit_spatial_boundary(
  df = meuse,
  y1 = "zinc",
  y2 = "lead",
  latent_col = "dist",
  iterations = 10000
)

cat("Acceptance rate:", round(fit$acceptance_rate * 100, 1), "%\n")
cat("Optimal c:", round(fit$optimal_c, 4), "\n")
cat("Trace head c:\n")
print(head(fit$fit$c, 20))
cat("Trace tail c:\n")
print(tail(fit$fit$c, 20))
