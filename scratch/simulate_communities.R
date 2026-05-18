library(ggplot2)
library(patchwork)

#' Evolutionary Particle Swarm for Spatial Communities
#' 
#' @param n_particles Total number of particles (divided equally into Red, Green, Blue)
#' @param clustering_strength How strongly particles of the SAME color attract each other
#' @param interaction_strength How strongly the cross-color rules (Red-Green affinity, Blue-Green repulsion) are applied
#' @param iterations Number of evolutionary steps
simulate_communities <- function(n_particles = 900, clustering_strength = 1.0, interaction_strength = 2.0, iterations = 30) {
  # Initialize particles
  types <- rep(c("Red", "Green", "Blue"), each = n_particles / 3)
  x <- runif(n_particles, 0, 100)
  y <- runif(n_particles, 0, 100)
  
  for (iter in 1:iterations) {
    # Compute pairwise distances (fast vectorized in R)
    dx <- outer(x, x, "-")
    dy <- outer(y, y, "-")
    dist_sq <- dx^2 + dy^2
    
    # We only care about local neighbors (radius of 15 units)
    # Convert distance to an attraction weight (closer = stronger pull)
    weight <- exp(-dist_sq / (2 * 10^2))
    diag(weight) <- 0 # Ignore self
    
    # Calculate Center of Mass pulls
    # Particles are pulled toward the weighted average position of their neighbors
    pull_x <- (weight %*% x) / (rowSums(weight) + 1e-5)
    pull_y <- (weight %*% y) / (rowSums(weight) + 1e-5)
    
    # Initialize movement vectors
    move_x <- numeric(n_particles)
    move_y <- numeric(n_particles)
    
    # 1. CLUSTERING (Self-Affinity)
    # Every particle wants to move toward its own color
    for (color in c("Red", "Green", "Blue")) {
      idx <- which(types == color)
      w_self <- weight[idx, idx]
      mx <- (w_self %*% x[idx]) / (rowSums(w_self) + 1e-5) - x[idx]
      my <- (w_self %*% y[idx]) / (rowSums(w_self) + 1e-5) - y[idx]
      move_x[idx] <- move_x[idx] + mx * clustering_strength
      move_y[idx] <- move_y[idx] + my * clustering_strength
    }
    
    # 2. EVOLUTIONARY RULES (Interaction)
    idx_R <- which(types == "Red")
    idx_G <- which(types == "Green")
    idx_B <- which(types == "Blue")
    
    # Red has Affinity for Green (Red moves TOWARD Green)
    w_RG <- weight[idx_R, idx_G]
    mx_RG <- (w_RG %*% x[idx_G]) / (rowSums(w_RG) + 1e-5) - x[idx_R]
    my_RG <- (w_RG %*% y[idx_G]) / (rowSums(w_RG) + 1e-5) - y[idx_R]
    move_x[idx_R] <- move_x[idx_R] + mx_RG * interaction_strength
    move_y[idx_R] <- move_y[idx_R] + my_RG * interaction_strength
    
    # Green has Affinity for Red (Green moves TOWARD Red)
    w_GR <- weight[idx_G, idx_R]
    mx_GR <- (w_GR %*% x[idx_R]) / (rowSums(w_GR) + 1e-5) - x[idx_G]
    my_GR <- (w_GR %*% y[idx_R]) / (rowSums(w_GR) + 1e-5) - y[idx_G]
    move_x[idx_G] <- move_x[idx_G] + mx_GR * interaction_strength
    move_y[idx_G] <- move_y[idx_G] + my_GR * interaction_strength
    
    # Blue has Repulsion for Green (Blue moves AWAY FROM Green)
    w_BG <- weight[idx_B, idx_G]
    mx_BG <- (w_BG %*% x[idx_G]) / (rowSums(w_BG) + 1e-5) - x[idx_B]
    my_BG <- (w_BG %*% y[idx_G]) / (rowSums(w_BG) + 1e-5) - y[idx_B]
    move_x[idx_B] <- move_x[idx_B] - mx_BG * interaction_strength
    move_y[idx_B] <- move_y[idx_B] - my_BG * interaction_strength
    
    # Update positions (with a small random walk to prevent getting stuck)
    x <- x + move_x * 0.2 + rnorm(n_particles, 0, 1.0)
    y <- y + move_y * 0.2 + rnorm(n_particles, 0, 1.0)
    
    # Clamp to boundaries
    x <- pmax(0, pmin(100, x))
    y <- pmax(0, pmin(100, y))
  }
  
  df <- data.frame(x = x, y = y, Color = types)
  return(df)
}

# Generate 4 different examples
set.seed(42)

message("Simulating Example 1...")
df1 <- simulate_communities(n_particles = 900, clustering_strength = 1.0, interaction_strength = 2.0, iterations = 30)
p1 <- ggplot(df1, aes(x, y, color=Color)) + geom_point(size=2, alpha=0.8) + 
  scale_color_manual(values=c("Red"="red", "Green"="#00FF00", "Blue"="dodgerblue")) +
  theme_void() + ggtitle("High Clustering, Strong Rules") + theme(legend.position="none", plot.title=element_text(hjust=0.5))

message("Simulating Example 2...")
df2 <- simulate_communities(n_particles = 900, clustering_strength = 0.5, interaction_strength = 3.0, iterations = 30)
p2 <- ggplot(df2, aes(x, y, color=Color)) + geom_point(size=2, alpha=0.8) + 
  scale_color_manual(values=c("Red"="red", "Green"="#00FF00", "Blue"="dodgerblue")) +
  theme_void() + ggtitle("Low Clustering, Aggressive Rules") + theme(legend.position="none", plot.title=element_text(hjust=0.5))

message("Simulating Example 3...")
df3 <- simulate_communities(n_particles = 900, clustering_strength = 2.0, interaction_strength = 0.5, iterations = 30)
p3 <- ggplot(df3, aes(x, y, color=Color)) + geom_point(size=2, alpha=0.8) + 
  scale_color_manual(values=c("Red"="red", "Green"="#00FF00", "Blue"="dodgerblue")) +
  theme_void() + ggtitle("Massive Clusters, Weak Rules") + theme(legend.position="none", plot.title=element_text(hjust=0.5))

message("Simulating Example 4...")
df4 <- simulate_communities(n_particles = 900, clustering_strength = 1.2, interaction_strength = 1.5, iterations = 30)
p4 <- ggplot(df4, aes(x, y, color=Color)) + geom_point(size=2, alpha=0.8) + 
  scale_color_manual(values=c("Red"="red", "Green"="#00FF00", "Blue"="dodgerblue")) +
  theme_void() + ggtitle("Balanced Evolution") + theme(legend.position="none", plot.title=element_text(hjust=0.5))

# Combine plots
final_plot <- (p1 | p2) / (p3 | p4) + plot_annotation(
  title = "Evolutionary Particle Swarm: Spatial Communities",
  subtitle = "Red=Affinity for Green | Blue=Repulsion from Green | Red/Blue=Independent",
  theme = theme(plot.title = element_text(size=18, face="bold"))
)

ggsave("simulated_communities.png", final_plot, width=10, height=10, bg="white")
message("Saved to simulated_communities.png")
