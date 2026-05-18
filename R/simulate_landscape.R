#' Simulate Alien Landscape
#' 
#' Generates a 2D grid representing continuous intensities of three microorganism colonies: Red, Green, and Blue.
#' 
#' Rules:
#' 1. Green and Blue boundaries are strictly disjoint (Repulsion).
#' 2. Red shares boundaries with Green in some regions (Affinity).
#' 3. Red shares boundaries with Blue in other regions (Co-habitation).
#' 
#' @param n_grid Integer, the size of the N x N grid.
#' @param noise Numeric, the standard deviation of the Gaussian noise added to the intensities.
#' @return A data.frame with columns x, y, Red, Green, and Blue.
#' @export
simulate_landscape <- function(n_grid = 100, noise = 0.5) {
  # Create the coordinate grid
  grid <- expand.grid(x = seq(0, 10, length.out = n_grid),
                      y = seq(0, 10, length.out = n_grid))
  
  # Initialize continuous background intensities (base levels)
  grid$Green <- 0
  grid$Blue <- 0
  grid$Red <- 0
  
  # --- 1. Define Green Colony (e.g., a circle on the left) ---
  # Center (3, 5), radius 2
  dist_to_green <- sqrt((grid$x - 3)^2 + (grid$y - 5)^2)
  in_green <- dist_to_green <= 2
  
  # Green surface: smooth dome inside the boundary
  grid$Green[in_green] <- 5 * (1 - dist_to_green[in_green]/2)
  
  # --- 2. Define Blue Colony (e.g., a circle on the right) ---
  # Center (7, 5), radius 2. Strictly disjoint from Green.
  dist_to_blue <- sqrt((grid$x - 7)^2 + (grid$y - 5)^2)
  in_blue <- dist_to_blue <= 2
  
  # Blue surface: smooth dome inside the boundary
  grid$Blue[in_blue] <- 5 * (1 - dist_to_blue[in_blue]/2)
  
  # --- 3. Define Red Colony ---
  # Red bridges both. It shares the Green boundary on the left, 
  # and the Blue boundary on the right, plus an area connecting them.
  # Let's make Red an ellipse that encapsulates both, so it exactly 
  # shares their outer boundaries but occupies the space between them too.
  
  # To make Red share exact boundary segments:
  # Red exists if it is inside Green OR inside Blue OR in a connecting bridge.
  # Let's create a rectangular bridge connecting the two circles.
  in_bridge <- (grid$x >= 3 & grid$x <= 7) & (grid$y >= 4 & grid$y <= 6)
  
  in_red <- in_green | in_blue | in_bridge
  
  # Red surface: a constant elevated plateau with noise
  grid$Red[in_red] <- 4
  
  # Add spatial noise to simulate real satellite sensor data
  grid$Green <- grid$Green + rnorm(nrow(grid), mean = 0, sd = noise)
  grid$Blue <- grid$Blue + rnorm(nrow(grid), mean = 0, sd = noise)
  grid$Red <- grid$Red + rnorm(nrow(grid), mean = 0, sd = noise)
  
  # Ensure no negative intensities (sensors don't read negative light)
  grid$Green <- pmax(0, grid$Green)
  grid$Blue <- pmax(0, grid$Blue)
  grid$Red <- pmax(0, grid$Red)
  
  return(grid)
}

#' Simulate Cortical Histological Tissue (Microglia, Astrocytes, PNNs)
#' 
#' Generates a 2D grid simulating biological spatial distributions of Microglia (Iba1),
#' Astrocytes (GFAP), and Perineuronal Networks (PNNs, e.g. WFA) under reactive neuroinflammation.
#' 
#' Biological Rules:
#' 1. Reactive Astrocytes (GFAP) form a dense inflammatory hotspot in the center (gliosis).
#' 2. Activated Microglia (Iba1) co-localize with the astrocytes (Affinity) in the gliosis zone.
#' 3. Perineuronal Networks (PNNs) wrap around neurons on the periphery but are actively
#'    degraded or excluded from the central inflammatory gliosis zone (Repulsion).
#' 
#' @param n_grid Integer, the size of the N x N grid.
#' @param noise Numeric, standard deviation of Gaussian noise.
#' @return A data.frame with columns x, y, Astrocytes, Microglia, PNNs, and LatentInflame.
#' @export
simulate_cortical_tissue <- function(n_grid = 100, noise = 0.3) {
  # Grid coordinates from 0 to 10
  grid <- expand.grid(x = seq(0, 10, length.out = n_grid),
                      y = seq(0, 10, length.out = n_grid))
  
  # 1. Define Central Inflammatory Latent Field (Gliosis gradient)
  # High in the center (5, 5)
  grid$LatentInflame <- 5 - sqrt((grid$x - 5)^2 + (grid$y - 5)^2)
  
  # Astrocytes (GFAP) trigger at a threshold of c = 2.0
  grid$Astrocytes <- 0.5 + 3.0 * pmax(0, grid$LatentInflame - 2.0)
  
  # Microglia (Iba1) share this boundary (Affinity) but have higher slope
  grid$Microglia <- 0.2 + 4.5 * pmax(0, grid$LatentInflame - 2.0)
  
  # 2. Define Perineuronal Networks (PNNs)
  # PNNs wrap around parvalbumin neurons located in the corners/periphery:
  # Neuron 1 at (2.5, 2.5), Neuron 2 at (7.5, 7.5)
  dist_n1 <- sqrt((grid$x - 2.5)^2 + (grid$y - 2.5)^2)
  dist_n2 <- sqrt((grid$x - 7.5)^2 + (grid$y - 7.5)^2)
  
  # PNNs exist in ring-like shells around these centers
  in_ring1 <- dist_n1 >= 0.8 & dist_n1 <= 1.8
  in_ring2 <- dist_n2 >= 0.8 & dist_n2 <= 1.8
  
  grid$PNNs <- 0.1
  grid$PNNs[in_ring1] <- 4.0 * (1.0 - abs(dist_n1[in_ring1] - 1.3)/0.5)
  grid$PNNs[in_ring2] <- 4.0 * (1.0 - abs(dist_n2[in_ring2] - 1.3)/0.5)
  
  # Inflammatory microglia degrade PNNs in the center:
  # PNNs are suppressed where inflammation is high
  grid$PNNs <- grid$PNNs * (1.0 - pmin(1.0, pmax(0, grid$LatentInflame - 1.5)/2.0))
  
  # Add noise
  grid$Astrocytes <- pmax(0, grid$Astrocytes + rnorm(nrow(grid), 0, noise))
  grid$Microglia <- pmax(0, grid$Microglia + rnorm(nrow(grid), 0, noise))
  grid$PNNs <- pmax(0, grid$PNNs + rnorm(nrow(grid), 0, noise))
  
  return(grid)
}

#' Calculate shortest distance to a line segment
#' @keywords internal
distance_to_segment <- function(px, py, ax, ay, bx, by) {
  ab_x <- bx - ax
  ab_y <- by - ay
  
  ap_x <- px - ax
  ap_y <- py - ay
  
  ab_len_sq <- ab_x^2 + ab_y^2
  if (ab_len_sq == 0) {
    return(sqrt(ap_x^2 + ap_y^2))
  }
  
  t <- (ap_x * ab_x + ap_y * ab_y) / ab_len_sq
  t <- pmax(0, pmin(1, t))
  
  cx <- ax + t * ab_x
  cy <- ay + t * ab_y
  
  sqrt((px - cx)^2 + (py - cy)^2)
}

#' Simulate Organic Branching Neuron Tissue
#' 
#' Generates a 2D grid simulating a highly realistic, organic, branching neuron
#' with a central soma and multiple bifurcating (forked) dendritic trees.
#' 
#' @param n_grid Integer, the size of the N x N grid.
#' @param noise Numeric, standard deviation of Gaussian noise.
#' @return A data.frame with columns x, y, Neuron, Microglia, PNNs, and LatentDendrite.
#' @export
simulate_dendritic_neuron <- function(n_grid = 100, noise = 0.2) {
  grid <- expand.grid(x = seq(0, 10, length.out = n_grid),
                      y = seq(0, 10, length.out = n_grid))
  
  # Define dendritic segments (Start X, Start Y, End X, End Y)
  segments <- list(
    # Primary Dendrite Trunk 1 (Up)
    list(5.0, 5.0, 5.0, 7.5),
    # Fork 1 (Left & Right)
    list(5.0, 7.5, 3.0, 9.0),
    list(5.0, 7.5, 7.0, 9.0),
    
    # Primary Dendrite Trunk 2 (Left Down)
    list(5.0, 5.0, 2.5, 3.5),
    # Fork 2
    list(2.5, 3.5, 0.8, 3.5),
    list(2.5, 3.5, 2.0, 1.2),
    
    # Primary Dendrite Trunk 3 (Right Down)
    list(5.0, 5.0, 7.5, 3.5),
    # Fork 3
    list(7.5, 3.5, 9.2, 3.5),
    list(7.5, 3.5, 8.0, 1.2)
  )
  
  # Calculate shortest distance to any segment for every grid point
  min_dist <- rep(Inf, nrow(grid))
  
  for (seg in segments) {
    dist_seg <- distance_to_segment(grid$x, grid$y, seg[[1]], seg[[2]], seg[[3]], seg[[4]])
    min_dist <- pmin(min_dist, dist_seg)
  }
  
  # Create a continuous Latent Dendritic Field (high near the dendrites, low far away)
  # Using a decaying exponential / Gaussian-like profile
  grid$LatentDendrite <- 10.0 * exp(-min_dist^2 / (2 * 0.8^2))
  
  # 1. Neuron Channel (e.g. MAP2)
  # Active near the dendrites (LatentDendrite > 2.0)
  grid$Neuron <- 0.2 + 4.5 * pmax(0, grid$LatentDendrite - 2.0)
  
  # 2. Microglia Channel (Iba1) showing Affinity
  # Microglia are active along the dendritic tree, wrapping synapses/branches (Affinity)
  grid$Microglia <- 0.1 + 3.0 * pmax(0, grid$LatentDendrite - 2.0)
  
  # 3. Perineuronal Networks (WFA) showing Repulsion / Spatial Exclusion
  # PNNs wrap primarily the soma (center (5, 5)) and proximal trunks, 
  # but are excluded from the active distal dendritic branches.
  # Let's say PNNs are active inside a sphere around the soma:
  dist_soma <- sqrt((grid$x - 5)^2 + (grid$y - 5)^2)
  pnn_profile <- 4.0 * exp(-dist_soma^2 / (2 * 1.5^2))
  
  # In addition, they are excluded from the active microglia/synaptic forks:
  # So they show Repulsion to the active dendritic field
  grid$PNNs <- pnn_profile * (1.0 - pmin(1.0, grid$Neuron/4.0))
  
  # Add noise
  grid$Neuron <- pmax(0, grid$Neuron + rnorm(nrow(grid), 0, noise))
  grid$Microglia <- pmax(0, grid$Microglia + rnorm(nrow(grid), 0, noise))
  grid$PNNs <- pmax(0, grid$PNNs + rnorm(nrow(grid), 0, noise))
  
  return(grid)
}

