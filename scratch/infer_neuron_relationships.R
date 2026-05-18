devtools::load_all()

# 1. Simulate the dendritic neuron channels
message("Simulating organic branching neuron histology...")
landscape <- simulate_dendritic_neuron(n_grid = 100, noise = 0.15)

# 2. Pack the stains into a list
stains <- list(
  Neuron = landscape$Neuron,
  Microglia = landscape$Microglia,
  PNNs = landscape$PNNs
)

# 3. Discover relationships automatically using our R-Rust Bayesian inference engine!
message("\nRunning Automatic Bayesian Relationship Discovery on organic structures...")
results <- discover_spatial_relationships(
  stains = stains,
  latent_field = landscape$LatentDendrite,
  iterations = 2000,   # 2000 iterations per chain
  bf_threshold = 10    # Strong evidence threshold
)

# 4. Print out the discovered relationships
cat("\n======================================================\n")
cat("          BAYESIAN INFERENCE RESULTS                  \n")
cat("======================================================\n")
print(results$relationships)
cat("======================================================\n")
