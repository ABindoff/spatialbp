#' Compute Latent Distance Field
#' 
#' Calculates the continuous nearest-neighbor distance field for a grid of points
#' against a set of target "source" locations. This converts an organic point process 
#' into a continuous spatial latent variable suitable for spatialbp.
#' 
#' @param grid_x Numeric vector of X coordinates for the evaluation grid.
#' @param grid_y Numeric vector of Y coordinates for the evaluation grid.
#' @param source_x Numeric vector of X coordinates for the source patches/communities.
#' @param source_y Numeric vector of Y coordinates for the source patches/communities.
#' @return A numeric vector of the same length as grid_x, representing the Euclidean distance to the nearest source point.
#' @export
compute_latent_distance <- function(grid_x, grid_y, source_x, source_y) {
  if (length(grid_x) != length(grid_y)) stop("grid_x and grid_y must have the same length.")
  if (length(source_x) != length(source_y)) stop("source_x and source_y must have the same length.")
  if (length(source_x) == 0) return(rep(Inf, length(grid_x)))
  
  # Ensure FNN is available (or we could use standard dist, but FNN is faster for large grids)
  if (!requireNamespace("FNN", quietly = TRUE)) {
    # Fallback to base R
    distances <- numeric(length(grid_x))
    for (i in seq_along(grid_x)) {
      dx <- source_x - grid_x[i]
      dy <- source_y - grid_y[i]
      distances[i] <- min(sqrt(dx^2 + dy^2))
    }
    return(distances)
  }
  
  grid_pts <- cbind(grid_x, grid_y)
  source_pts <- cbind(source_x, source_y)
  
  nn <- FNN::get.knnx(data = source_pts, query = grid_pts, k = 1)
  return(as.numeric(nn$nn.dist))
}

#' Compute Latent Density Field
#' 
#' Calculates a continuous 2D Kernel Density Estimate (KDE) field for a set of points.
#' 
#' @param grid_x Numeric vector of X coordinates for the full output grid (e.g., from expand.grid).
#' @param grid_y Numeric vector of Y coordinates for the full output grid.
#' @param source_x Numeric vector of X coordinates for the source points.
#' @param source_y Numeric vector of Y coordinates for the source points.
#' @param res Integer, the resolution of the internal KDE grid (default 100).
#' @return A numeric vector of densities perfectly mapped to the (grid_x, grid_y) pairs.
#' @export
compute_latent_density <- function(grid_x, grid_y, source_x, source_y, res = 100) {
  if (length(grid_x) != length(grid_y)) stop("grid_x and grid_y must have the same length.")
  if (length(source_x) < 2) return(rep(0, length(grid_x)))
  
  xlim <- range(grid_x, na.rm = TRUE)
  ylim <- range(grid_y, na.rm = TRUE)
  
  # Add small buffer to limits to prevent edge artifacts
  xlim <- xlim + c(-1, 1) * diff(xlim) * 0.05
  ylim <- ylim + c(-1, 1) * diff(ylim) * 0.05
  
  # Compute KDE
  kde <- MASS::kde2d(source_x, source_y, n = res, lims = c(xlim, ylim))
  
  # We need to map the internal `kde` grid back to the arbitrary `grid_x, grid_y` vectors.
  # We can use bilinear interpolation via `fields::interp.surface` if available.
  if (requireNamespace("fields", quietly = TRUE)) {
    kde_obj <- list(x = kde$x, y = kde$y, z = kde$z)
    mapped_densities <- fields::interp.surface(kde_obj, cbind(grid_x, grid_y))
    mapped_densities[is.na(mapped_densities)] <- 0
    return(mapped_densities * 1000) # Scale for numerical stability in MCMC
  } else {
    # Slower fallback mapping
    mapped_densities <- numeric(length(grid_x))
    for (i in seq_along(grid_x)) {
      ix <- which.min(abs(kde$x - grid_x[i]))
      iy <- which.min(abs(kde$y - grid_y[i]))
      mapped_densities[i] <- kde$z[ix, iy]
    }
    return(mapped_densities * 1000)
  }
}
