#' Spatial Hinge Function (Continuous Level Set)
#' 
#' The mathematical core of 2D piecewise regression. Applies a threshold `c` 
#' to a latent spatial field `z`. It enforces spatial continuity (the surfaces
#' stitch together exactly at the boundary contour).
#' 
#' @param z Numeric vector representing the latent spatial propensity field.
#' @param c Numeric scalar, the breakpoint/threshold contour.
#' @return A numeric vector of the hinged basis $(z - c)_+$
#' @export
spatial_hinge <- function(z, c) {
  # The max(0, z - c) operation. Turns "on" only when z exceeds the threshold.
  pmax(0, z - c)
}

#' Spatial Mask (Discontinuous Level Set)
#' 
#' A stricter boundary where the surface is allowed to abruptly jump
#' (discontinuous) rather than smoothly hinge.
#' 
#' @param z Numeric vector representing the latent spatial propensity field.
#' @param c Numeric scalar, the breakpoint/threshold contour.
#' @return An integer vector (0 or 1) spatial mask.
#' @export
level_set_mask <- function(z, c) {
  as.integer(z > c)
}

#' Generate Radial Basis Functions (Thin-plate splines)
#' 
#' Generates the underlying smooth spatial surface Z(x,y) from a set of knots.
#' In a Bayesian framework, the weights of these basis functions are estimated
#' to warp the boundary contours around the protein puncta.
#' 
#' @param x Numeric vector of x coordinates
#' @param y Numeric vector of y coordinates
#' @param knots_x Numeric vector of knot x coordinates
#' @param knots_y Numeric vector of knot y coordinates
#' @return A matrix of radial basis evaluations
#' @export
construct_rbf <- function(x, y, knots_x, knots_y) {
  n <- length(x)
  k <- length(knots_x)
  
  # Distance matrix between data points and knots
  basis_matrix <- matrix(0, nrow = n, ncol = k)
  
  for (i in seq_len(k)) {
    # Euclidean distance
    r <- sqrt((x - knots_x[i])^2 + (y - knots_y[i])^2)
    
    # Thin-plate spline radial basis function: r^2 * log(r)
    # Handle the limit as r -> 0 where r^2 * log(r) = 0
    tps <- r^2 * log(r)
    tps[r == 0] <- 0 
    
    basis_matrix[, i] <- tps
  }
  
  return(basis_matrix)
}
