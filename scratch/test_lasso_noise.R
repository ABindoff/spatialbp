library(spatialbp)

data(meuse, package = "sp")

set.seed(42)
meuse$noise_channel <- rnorm(nrow(meuse), mean = 0, sd = 1)

# Fit without L1 regularization (lambda = 0)
message("Fitting WITH L1 Regularization (lambda = 0.0)...")
fit_no_reg <- fit_spatial_boundary(
  df = meuse,
  y1 = "zinc",
  y2 = "noise_channel",
  latent_col = "dist",
  iterations = 5000,
  lambda = 0.0
)

# Fit WITH large L1 regularization (lambda = 50.0)
message("\nFitting WITH L1 Regularization (lambda = 50.0)...")
fit_reg <- fit_spatial_boundary(
  df = meuse,
  y1 = "zinc",
  y2 = "noise_channel",
  latent_col = "dist",
  iterations = 5000,
  lambda = c(0.0, 50.0) # Penalty is 50!
)

cat("\n--- Standardized Beta2 (Noise Channel) Trace Summary ---\n")
cat("Without Regularization:\n")
print(summary(fit_no_reg$fit$beta2[2000:5000]))

cat("\nWith Regularization (lambda = 50.0):\n")
print(summary(fit_reg$fit$beta2[2000:5000]))
