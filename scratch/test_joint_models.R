devtools::load_all()

# 1. Simulate the landscape
message("Simulating landscape...")
landscape <- simulate_landscape(n_grid = 30, noise = 0.2)

# 2. Define the continuous latent spatial field
center_x <- 5
center_y <- 5
latent_field <- 10 - sqrt((landscape$x - center_x)^2 + (landscape$y - center_y)^2)

# 3. Fit and test the joint models
message("Running Multivariate Spatial Affinity and Repulsion Tests...")
results <- test_spatialbp(
  red = landscape$Red,
  green = landscape$Green,
  blue = landscape$Blue,
  latent_field = latent_field,
  iterations = 2000
)

# 4. Print results
cat("\n--- Joint Hypothesis Testing Results ---\n")
cat("Independent Null Model Log-Likelihood:", round(results$ll_null, 2), "\n")
cat("Joint Affinity Model (Red & Green) Log-Likelihood:", round(results$ll_affinity, 2), "\n")
cat("Joint Repulsion Model (Green & Blue) Log-Likelihood:", round(results$ll_repulsion, 2), "\n")
cat("\n--- Bayes Factors ---\n")
cat("Bayes Factor for Affinity (Red & Green):", results$bayes_factor_affinity, "\n")
cat("Bayes Factor for Repulsion (Green & Blue):", results$bayes_factor_repulsion, "\n")
