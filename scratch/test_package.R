library(spatialbp)

# 1. Simulate the landscape
message("Simulating landscape...")
landscape <- simulate_landscape(n_grid = 30, noise = 0.2) # smaller grid for fast test

# 2. Define a simple latent spatial field
# Let's say the latent field is a distance field from the center of the Red colony bridge
center_x <- 5
center_y <- 5
latent_field <- 10 - sqrt((landscape$x - center_x)^2 + (landscape$y - center_y)^2)

# 3. Fit the spatial piecewise model for Red using the latent field
message("Fitting spatialbp model using Rust MCMC engine...")
fit <- fit_spatialbp(
  observed = landscape$Red,
  baseline = rep(0, nrow(landscape)),
  latent = latent_field,
  iterations = 1000
)

# 4. Print results
message("\n--- Fit Summary ---")
cat("Acceptance rate:", fit$acceptance_rate, "\n")
cat("Posterior median of beta:", median(fit$traces$beta), "\n")
cat("Posterior median of cp (changepoint):", median(fit$traces$cp), "\n")
cat("Log-likelihood range:", range(fit$log_likelihood), "\n")
