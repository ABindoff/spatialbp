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
