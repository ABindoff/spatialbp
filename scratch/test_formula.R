devtools::load_all()

# Create dummy dataset
df <- data.frame(
  Green = rnorm(100, 10, 2),
  Red = rnorm(100, 5, 1),
  radial_z = runif(100, 0, 100),
  tissue_thickness = rnorm(100, 20, 5),
  oxygen = runif(100, 0, 1),
  inflammation = rbinom(100, 1, 0.5)
)

message("\n--- Testing 1D Formula with Baseline ---")
fit1 <- spatialbp(
  Green ~ tissue_thickness + hinge(latent = radial_z),
  data = df
)

message("\n--- Testing Conditional Modifiers ---")
fit2 <- spatialbp(
  Green ~ 1 + hinge(latent = radial_z, mod_beta = ~ oxygen, mod_c = ~ inflammation, smooth = TRUE),
  data = df
)

message("Summary of Fit 2:")
print(str(fit2$traces))

# Create a simple diagnostic plot
png("spatialbp_formula_diagnostics.png", width=800, height=400)
par(mfrow=c(1,2))
plot(fit2$traces$beta, type="l", main="Trace: Beta Intercept", ylab="Beta")
plot(fit2$traces$cp, type="l", main="Trace: Boundary (c) Intercept", ylab="c")
dev.off()

message("\n--- Testing Multivariate + Offset ---")
fit3 <- spatialbp(
  cbind(Green, Red) ~ offset(tissue_thickness) + hinge(latent = radial_z, mod_beta = ~ oxygen),
  data = df
)
