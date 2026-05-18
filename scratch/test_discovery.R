devtools::load_all()

# 1. Simulate the landscape
message("Simulating landscape...")
landscape <- simulate_landscape(n_grid = 30, noise = 0.2)

# 2. Define the X coordinate as the latent field (partitioning left vs right)
latent_field <- landscape$x

# 3. Pack the stains into a list
stains <- list(
  Red = landscape$Red,
  Green = landscape$Green,
  Blue = landscape$Blue
)

# 4. Discover relationships automatically!
message("Running Automatic Spatial Relationship Discovery...")
results <- discover_spatial_relationships(
  stains = stains,
  latent_field = latent_field,
  iterations = 1000,   # 1000 iterations per fit for speed
  bf_threshold = 3.0   # Positive evidence threshold
)

# 5. Print out the discovered relationships
cat("\n--- Discovered Pairwise Spatial Relationships ---\n")
print(results$relationships)
