#' Fit a 2D Spatial Piecewise Regression Model
#' 
#' Wraps the Rust MCMC engine and structures the output to perfectly mimic
#' a `smoothbp` object, inheriting all posterior diagnostic tools.
#' 
#' @param observed Numeric vector of observed intensities
#' @param baseline Numeric vector of baseline values
#' @param latent Numeric vector of latent spatial field values
#' @param iterations Integer, number of MCMC iterations
#' @return An object of class \code{c("spatialbp", "smoothbp")}
#' @export
fit_spatialbp <- function(observed, baseline, latent, iterations = 5000) {
  
  # 1. Call the Rust MCMC Engine
  raw_fit <- run_spatial_mcmc(
    observed = observed,
    baseline = baseline,
    latent = latent,
    iterations = iterations,
    init_beta = 1.0, init_c = 0.0,
    sd_beta = 0.1, sd_c = 0.1
  )
  
  # 2. Standardize Parameter Names for smoothbp 
  # We rename the 2D threshold 'c' to 'cp' (changepoint) so that smoothbp's
  # summary() and loo() methods can seamlessly discover it.
  formatted_fit <- list(
    traces = list(
      beta = raw_fit$beta,
      cp = raw_fit$c 
    ),
    log_likelihood = raw_fit$log_likelihood,
    acceptance_rate = raw_fit$acceptance_rate,
    
    # We must store the spatial context because smoothbp doesn't expect 2D data
    spatial_context = list(
      observed = observed,
      latent = latent
    )
  )
  
  # 3. Apply the Object-Oriented Inheritance
  # This makes R check for spatialbp functions first, then fallback to smoothbp
  class(formatted_fit) <- c("spatialbp", "smoothbp")
  
  return(formatted_fit)
}

#' Confidence Intervals for Spatial Boundaries
#' 
#' Computes the Highest Density Interval (HDI) or standard quantiles for the estimated 
#' boundary parameter `c` and other MCMC traces in a spatialbp model.
#' 
#' @param object A spatialbp model object.
#' @param parm Character vector of parameters to calculate intervals for (default: "cp").
#' @param level Numeric, confidence level (default: 0.95).
#' @param ... Unused.
#' @importFrom stats quantile
#' @export
confint.spatialbp <- function(object, parm = "cp", level = 0.95, ...) {
  if (is.null(object$traces)) {
    stop("No MCMC traces found in spatialbp object.")
  }
  
  alpha <- 1 - level
  lower_q <- alpha / 2
  upper_q <- 1 - (alpha / 2)
  
  if (missing(parm) || is.null(parm)) {
    parm <- "cp"
  }
  
  results <- matrix(NA, nrow = length(parm), ncol = 2)
  rownames(results) <- parm
  colnames(results) <- c(paste0(lower_q * 100, "%"), paste0(upper_q * 100, "%"))
  
  for (i in seq_along(parm)) {
    p <- parm[i]
    if (p %in% names(object$traces)) {
      results[i, ] <- stats::quantile(object$traces[[p]], probs = c(lower_q, upper_q), na.rm = TRUE)
    } else {
      warning("Parameter '", p, "' not found in MCMC traces.")
    }
  }
  
  return(results)
}

#' Plot a 2D Spatial Piecewise Regression
#' 
#' Intercepts the generic plot() call to prevent smoothbp from trying to draw
#' a 1D line. Instead, this renders a 2D topographical map of the boundaries.
#' 
#' @param z_slice Numeric (optional), specific Z-slice to display for 3D/4D data.
#' @param t_slice Numeric (optional), specific time-point to display for 4D data.
#' @param facet_by Character (optional), column name to facet by.
#' @param ... Unused.
#' @importFrom ggplot2 ggplot aes geom_raster geom_point geom_tile scale_fill_viridis_c scale_color_viridis_c labs theme_minimal facet_wrap
#' @export
plot.spatialbp <- function(x, z_slice = NULL, t_slice = NULL, facet_by = NULL, ...) {
  message("Rendering 2D/3D Spatial Boundary Map...")
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting.")
  }
  
  if (is.null(x$data)) {
    stop("Original data frame not found in model object. Plotting requires x$data.")
  }
  
  df <- x$data
  
  # Calculate post-warmup median boundary
  if (!is.null(x$warmup) && x$warmup > 0 && x$warmup < length(x$traces$cp)) {
    optimal_c <- median(x$traces$cp[(x$warmup + 1):length(x$traces$cp)])
  } else {
    optimal_c <- median(x$traces$cp)
  }
  
  latent_field <- x$spatial_context$latent
  # Use the first response channel for the background raster
  bg_channel <- colnames(x$spatial_context$observed)[1]
  
  # Calculate boundary spots (points close to estimated changepoint c)
  tolerance <- (max(latent_field) - min(latent_field)) * 0.025
  df$Boundary_Zone <- abs(latent_field - optimal_c) < tolerance
  
  # Check for 3D/4D coordinates
  has_z <- "z" %in% colnames(df)
  has_t <- "t" %in% colnames(df)
  
  if (has_t && !is.null(t_slice)) {
    unique_t <- unique(df$t)
    closest_t <- unique_t[which.min(abs(unique_t - t_slice))]
    df <- df[df$t == closest_t, ]
    message("Slicing time-point t = ", closest_t)
  }
  
  if (has_z && !is.null(z_slice) && is.null(facet_by)) {
    unique_z <- unique(df$z)
    closest_z <- unique_z[which.min(abs(unique_z - z_slice))]
    df <- df[df$z == closest_z, ]
    message("Slicing Z-plane z = ", closest_z)
  }
  
  if (has_z && is.null(z_slice) && is.null(facet_by)) {
    unique_z <- unique(df$z)
    median_z <- unique_z[order(unique_z)[length(unique_z) %/% 2 + 1]]
    df <- df[df$z == median_z, ]
    message("No Z-slice specified. Defaulting to central Z-plane z = ", round(median_z, 2))
  }
  
  unique_x <- length(unique(df$x))
  unique_y <- length(unique(df$y))
  is_grid <- (unique_x * unique_y <= nrow(df) * 1.3)
  
  boundary_color <- "cyan"
  
  if (is_grid) {
    step_x <- min(diff(sort(unique(df$x))))
    step_y <- min(diff(sort(unique(df$y))))
    if (is.infinite(step_x) || step_x == 0) step_x <- 1
    if (is.infinite(step_y) || step_y == 0) step_y <- 1
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = .data[[bg_channel]])) +
      ggplot2::scale_fill_viridis_c(option = "magma", name = bg_channel) +
      ggplot2::geom_tile(data = subset(df, Boundary_Zone), 
                         ggplot2::aes(x = x, y = y), 
                         width = step_x, height = step_y,
                         fill = boundary_color, color = boundary_color, alpha = 0.8)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(ggplot2::aes(color = .data[[bg_channel]]), size = 3, alpha = 0.8) +
      ggplot2::scale_color_viridis_c(option = "magma", name = bg_channel) +
      ggplot2::geom_point(data = subset(df, Boundary_Zone),
                          ggplot2::aes(x = x, y = y),
                          shape = 21, color = boundary_color, fill = NA, size = 5, stroke = 1.5)
  }
  
  if (!is.null(facet_by) && facet_by %in% colnames(df)) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_by)))
  }
  
  subtitle_text <- paste("Channel:", bg_channel, "\nBoundary c =", round(optimal_c, 3))
  if (has_z && is.null(facet_by)) subtitle_text <- paste(subtitle_text, "\nSlice: z =", round(df$z[1], 2))
  if (has_t && !is.null(t_slice)) subtitle_text <- paste(subtitle_text, ", t =", round(df$t[1], 2))
  
  p <- p + 
    ggplot2::labs(title = "2D Spatial Boundary Map", subtitle = subtitle_text) +
    ggplot2::theme_minimal()
  
  return(p)
}
