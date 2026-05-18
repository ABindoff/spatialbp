source("../R/simulate_landscape.R")
source("../R/spatial_basis.R")

cat("=== spatialbp End-to-End Pipeline Demo ===\n\n")

# 1. Generate the Simulated Data
set.seed(42)
cat("[1/4] Generating Landscape...\n")
landscape <- simulate_landscape(n_grid = 50, noise = 0.1)

# 2. Define Spatial Knots
# We lay down a grid of control points that will anchor our continuous surfaces.
knots_x <- seq(1, 9, length.out = 5)
knots_y <- seq(1, 9, length.out = 5)
knots <- expand.grid(x = knots_x, y = knots_y)

# 3. Construct the Latent Propensity Field (Z)
cat("[2/4] Constructing Thin-Plate Spline Radial Basis Matrix...\n")
basis_matrix <- construct_rbf(landscape$x, landscape$y, knots$x, knots$y)

# For this demo, we generate a synthetic latent field using arbitrary weights.
# In the full Rust MCMC, these weights are parameters that are sampled and optimized.
set.seed(123)
latent_weights <- rnorm(ncol(basis_matrix), 0, 1)
latent_field <- as.numeric(basis_matrix %*% latent_weights)

# 4. Prepare the baseline surface (U) 
# Assuming a flat baseline of 0 for simplicity
baseline <- rep(0, nrow(landscape))

cat("[3/4] Data structures initialized.\n")
cat(sprintf("      - Total Pixels: %d\n", nrow(landscape)))
cat(sprintf("      - Spatial Knots: %d\n", nrow(knots)))
cat(sprintf("      - Basis Matrix Size: %d x %d\n", nrow(basis_matrix), ncol(basis_matrix)))

# 5. The Rust MCMC Handoff
cat("\n[4/4] Ready for Rust MCMC Handoff.\n")
cat("      If the Rust DLL were loaded, we would now execute:\n")
cat("      --------------------------------------------------\n")
cat("      trace_green <- run_spatial_mcmc(\n")
cat("        observed = landscape$Green,\n")
cat("        baseline = baseline,\n")
cat("        latent = latent_field,\n")
cat("        iterations = 10000,\n")
cat("        init_beta = 1.0, init_c = 0.0,\n")
cat("        sd_beta = 0.1, sd_c = 0.1\n")
cat("      )\n")
cat("      --------------------------------------------------\n")
cat("Pipeline execution complete.\n")
