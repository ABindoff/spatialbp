#' High-Level Spatial Boundary Fitting
#' 
#' A biologist-friendly wrapper that automatically calibrates, runs, and fits
#' spatial boundary models without requiring manual tuning of MCMC proposal settings.
#' Exposes L1 regularization (Laplace priors) to filter out noise in multi-boundary datasets.
#' 
#' @param df A data.frame containing spatial coordinate columns (x, y) and protein intensities.
#' @param y1 Character, column name of the first protein channel.
#' @param y2 Character (optional), column name of the second protein channel for joint modeling.
#' @param latent_type Character, type of boundary geometry to construct: `"radial"` (circular boundary around center) or `"linear"` (linear boundary along X-axis). Default is `"radial"`.
#' @param center Numeric vector of length 2, the center coordinates for the radial boundary (optional, defaults to geometric center of the coordinates).
#' @param iterations Integer, number of MCMC iterations. Default is `3000`.
#' @param latent_col Character (optional), column name of a pre-calculated latent environmental gradient (e.g. `"dist"`, `"elev"`, or `"depth"`). If provided, `latent_type` is ignored.
#' @param lambda Numeric or numeric vector of length 2 (optional), Laplace regularization penalty parameter to shrink noise channels to 0. Default is `0.0`.
#' @return A `spatialbp_fit` object containing estimated boundaries and MCMC traces.
#' @export
fit_spatial_boundary <- function(df, y1, y2 = NULL, latent_type = "radial", center = NULL, iterations = 3000, latent_col = NULL, lambda = 0.0) {
  # 1. Coordinate checks
  if (!all(c("x", "y") %in% colnames(df))) {
    stop("Data frame must contain 'x' and 'y' coordinate columns.")
  }
  
  if (!y1 %in% colnames(df)) {
    stop(paste("Channel", y1, "not found in data frame."))
  }
  
  # 2. Define the continuous latent coordinate gradient
  if (!is.null(latent_col)) {
    if (!latent_col %in% colnames(df)) {
      stop(paste("Latent column", latent_col, "not found in data frame."))
    }
    latent_field <- df[[latent_col]]
    message("Using custom latent gradient column: ", latent_col)
  } else {
    if (latent_type == "radial") {
      if (is.null(center)) {
        center <- c(mean(range(df$x)), mean(range(df$y)))
      }
      # Z represents distance from the center (active close to center)
      latent_field <- 10 - sqrt((df$x - center[1])^2 + (df$y - center[2])^2)
      message("Constructed radial latent field around center (", round(center[1], 2), ", ", round(center[2], 2), ")")
    } else if (latent_type == "linear") {
      # Z represents horizontal position along the X-axis
      latent_field <- df$x
      message("Constructed linear latent field along X-axis.")
    } else {
      stop("Unknown latent_type. Must be 'radial' or 'linear'.")
    }
  }
  
  # 3. Handle gradient direction (Hinge-Up vs Hinge-Down) dynamically
  # If the response y1 decreases as the latent field increases, negate the latent field
  # to convert it into a standard Hinge-Up shape for the Rust MCMC backend!
  y1_orig <- as.double(df[[y1]])
  corr_y1_latent <- cor(y1_orig, latent_field, use = "complete.obs")
  negate_latent <- FALSE
  if (!is.na(corr_y1_latent) && corr_y1_latent < 0) {
    message("Detected Hinge-Down gradient direction. Negating latent coordinates internally for MCMC...")
    negate_latent <- TRUE
    latent_field_mcmc <- -latent_field
  } else {
    latent_field_mcmc <- latent_field
  }
  
  # 4. Adaptive Proposal Calibration
  message("Calibrating MCMC proposal step sizes automatically...")
  latent_range <- max(latent_field_mcmc, na.rm = TRUE) - min(latent_field_mcmc, na.rm = TRUE)
  if (latent_range == 0) latent_range <- 1
  
  # Mathematically scale step sizes by 1 / sqrt(N) to prevent the log-likelihood from being
  # extremely hyper-sensitive on large datasets.
  # We cap the scale factor at 15.0 so that massive raster grids don't shrink the step size
  # down to microscopic levels, guaranteeing fast convergence and a healthy "caterpillar" trace!
  n_samples <- nrow(df)
  scale_factor <- min(max(sqrt(n_samples), 5.0), 15.0)
  
  sd_beta <- (0.1 / latent_range) / scale_factor
  init_beta_val <- 1.0 / latent_range
  sd_c <- (latent_range * 0.15) / scale_factor
  
  # 5. Perform MCMC Estimation
  # Force all vectors to double and scale them to have mean = 0 and sd = 1
  mean_y1 <- mean(y1_orig, na.rm = TRUE)
  sd_y1 <- sd(y1_orig, na.rm = TRUE)
  if (sd_y1 == 0) sd_y1 <- 1
  y1_vec <- (y1_orig - mean_y1) / sd_y1
  
  latent_double <- as.double(latent_field_mcmc)
  
  # Mathematically rigorous baseline shifts: set baseline to the minimum scaled intensity
  # (so that the piecewise hinge represents growth above the minimum baseline!)
  y1_base_val <- min(y1_vec, na.rm = TRUE)
  baseline1 <- as.double(rep(y1_base_val, nrow(df)))
  
  if (is.null(y2)) {
    lambda_val <- if (is.null(lambda)) 0.0 else as.double(lambda[1])
    message("Fitting single-channel spatial boundary for: ", y1)
    fit <- run_spatial_mcmc(
      observed = y1_vec,
      baseline = baseline1,
      latent = latent_double,
      iterations = iterations,
      init_beta = init_beta_val,
      init_c = mean(latent_double, na.rm = TRUE),
      sd_beta = sd_beta,
      sd_c = sd_c,
      lambda = lambda_val
    )
    
    # Scale back MCMC trace and optimal parameter estimations
    if (negate_latent) {
      fit$c <- -fit$c
    }
    optimal_c <- median(fit$c)
    optimal_beta <- median(fit$beta) * sd_y1
    
    res <- list(
      y1 = y1,
      y2 = NULL,
      mode = "single",
      optimal_c = optimal_c,
      optimal_beta = optimal_beta,
      acceptance_rate = fit$acceptance_rate,
      latent_field = as.double(latent_field),
      df = df,
      fit = fit
    )
  } else {
    if (!y2 %in% colnames(df)) {
      stop(paste("Channel", y2, "not found in data frame."))
    }
    
    y2_orig <- as.double(df[[y2]])
    mean_y2 <- mean(y2_orig, na.rm = TRUE)
    sd_y2 <- sd(y2_orig, na.rm = TRUE)
    if (sd_y2 == 0) sd_y2 <- 1
    y2_vec <- (y2_orig - mean_y2) / sd_y2
    
    # Establish the base value for channel 2
    y2_base_val <- min(y2_vec, na.rm = TRUE)
    baseline2 <- as.double(rep(y2_base_val, nrow(df)))
    
    # Check if they co-localize or segregate to determine Affinity vs Repulsion automatically
    # By looking at their Pearson correlation coefficient of original scales:
    corr <- cor(y1_orig, y2_orig, use = "complete.obs")
    mode <- if (corr >= 0) "affinity" else "repulsion"
    
    message("Detected correlation of ", round(corr, 3), " between ", y1, " and ", y2)
    message("Automatically fitting Joint ", toupper(mode), " Model...")
    
    lambda_val1 <- if (is.null(lambda)) 0.0 else as.double(lambda[1])
    lambda_val2 <- if (is.null(lambda) || length(lambda) < 2) lambda_val1 else as.double(lambda[2])
    
    fit <- run_joint_mcmc(
      y1 = y1_vec,
      y2 = y2_vec,
      baseline1 = baseline1,
      baseline2 = baseline2,
      latent = latent_double,
      mode = mode,
      iterations = iterations,
      init_beta1 = init_beta_val, init_beta2 = init_beta_val,
      init_c = mean(latent_double, na.rm = TRUE),
      sd_beta1 = sd_beta, sd_beta2 = sd_beta,
      sd_c = sd_c,
      lambda1 = lambda_val1,
      lambda2 = lambda_val2
    )
    
    # Scale back MCMC trace and optimal parameter estimations
    if (negate_latent) {
      fit$c <- -fit$c
    }
    optimal_c <- median(fit$c)
    optimal_beta1 <- median(fit$beta1) * sd_y1
    optimal_beta2 <- median(fit$beta2) * sd_y2
    
    res <- list(
      y1 = y1,
      y2 = y2,
      mode = mode,
      optimal_c = optimal_c,
      optimal_beta1 = optimal_beta1,
      optimal_beta2 = optimal_beta2,
      acceptance_rate = fit$acceptance_rate,
      latent_field = as.double(latent_field),
      df = df,
      fit = fit
    )
  }
  
  class(res) <- "spatialbp_fit"
  message("Model fit successfully! Acceptance rate: ", round(res$acceptance_rate * 100, 1), "%")
  return(res)
}

#' Plot Spatial Boundary Fits
#' 
#' Automatically generates a premium, publication-grade overlay plot showing
#' the fitted continuous spatial boundary contour. Natively supports 2D, 3D, and 4D
#' datasets by automatically slicing or facetting across Z-planes or time-points.
#' Detects whether the input is a regular grid or sparse spatial point pattern,
#' rendering rasters or point-wise halos dynamically.
#' 
#' @param x A `spatialbp_fit` object.
#' @param z_slice Numeric (optional), specific Z-slice to display for 3D/4D data.
#' @param t_slice Numeric (optional), specific time-point to display for 4D data.
#' @param facet_by Character (optional), column name to facet by (e.g. `"z"` or `"t"`).
#' @param ... Unused.
#' @importFrom ggplot2 ggplot aes geom_raster geom_point geom_tile scale_fill_viridis_c scale_color_viridis_c labs theme_minimal facet_wrap geom_point
#' @export
plot.spatialbp_fit <- function(x, z_slice = NULL, t_slice = NULL, facet_by = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting.")
  }
  
  df <- x$df
  optimal_c <- x$optimal_c
  latent_field <- x$latent_field
  bg_channel <- x$y1
  
  # Calculate boundary spots (points close to estimated changepoint c)
  # Dynamic tolerance depending on grid spacing
  tolerance <- (max(latent_field) - min(latent_field)) * 0.025
  df$Boundary_Zone <- abs(latent_field - optimal_c) < tolerance
  
  # 1. Check for 3D/4D coordinates
  has_z <- "z" %in% colnames(df)
  has_t <- "t" %in% colnames(df)
  
  # Filter by time slice if 4D
  if (has_t && !is.null(t_slice)) {
    unique_t <- unique(df$t)
    closest_t <- unique_t[which.min(abs(unique_t - t_slice))]
    df <- df[df$t == closest_t, ]
    message("Slicing time-point t = ", closest_t)
  }
  
  # Filter by Z-slice if 3D/4D
  if (has_z && !is.null(z_slice) && is.null(facet_by)) {
    unique_z <- unique(df$z)
    closest_z <- unique_z[which.min(abs(unique_z - z_slice))]
    df <- df[df$z == closest_z, ]
    message("Slicing Z-plane z = ", closest_z)
  }
  
  # If 3D/4D but no slices or facets are specified, default to a smart central slice
  if (has_z && is.null(z_slice) && is.null(facet_by)) {
    unique_z <- unique(df$z)
    median_z <- unique_z[order(unique_z)[length(unique_z) %/% 2 + 1]]
    df <- df[df$z == median_z, ]
    message("No Z-slice specified. Defaulting to central Z-plane z = ", round(median_z, 2))
  }
  
  # 2. Smart Grid vs Point Pattern Detection
  # Regular grids have few unique coordinates relative to row count
  unique_x <- length(unique(df$x))
  unique_y <- length(unique(df$y))
  is_grid <- (unique_x * unique_y <= nrow(df) * 1.3)
  
  # Boundary color scheme (glowing green for segregation/repulsion, cyan for affinity)
  boundary_color <- if (x$mode == "repulsion") "chartreuse" else "cyan"
  
  if (is_grid) {
    # Dense regular grid layout
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = .data[[bg_channel]])) +
      ggplot2::scale_fill_viridis_c(option = "magma", name = bg_channel) +
      ggplot2::geom_tile(data = subset(df, Boundary_Zone), 
                         ggplot2::aes(x = x, y = y), 
                         fill = boundary_color, color = boundary_color, alpha = 0.8)
  } else {
    # Irregular sparse spatial point pattern
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(ggplot2::aes(color = .data[[bg_channel]]), size = 3, alpha = 0.8) +
      ggplot2::scale_color_viridis_c(option = "magma", name = bg_channel) +
      # Draw bright highlighted circular rings (halos) around points near the boundary
      ggplot2::geom_point(data = subset(df, Boundary_Zone),
                          ggplot2::aes(x = x, y = y),
                          shape = 21, color = boundary_color, fill = NA, size = 5, stroke = 1.5)
  }
  
  # Apply facetting if requested
  if (!is.null(facet_by) && facet_by %in% colnames(df)) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_by)))
  }
  
  # Dynamic labels
  subtitle_text <- paste("Channel:", x$y1, 
                         if(!is.null(x$y2)) paste("+", x$y2, "(Joint", toupper(x$mode), ")") else "",
                         "\nBoundary Changepoint c =", round(optimal_c, 3))
  
  if (has_z && is.null(facet_by)) {
    subtitle_text <- paste(subtitle_text, "\nSlice: z =", round(df$z[1], 2))
  }
  if (has_t && !is.null(t_slice)) {
    subtitle_text <- paste(subtitle_text, ", t =", round(df$t[1], 2))
  }
  
  p <- p + 
    ggplot2::labs(
      title = "Estimated Spatial Boundary Contour",
      subtitle = subtitle_text,
      x = "X Coordinate", y = "Y Coordinate"
    ) +
    ggplot2::theme_minimal()
  
  return(p)
}

#' Plot MCMC Diagnostic Traces
#' 
#' Generates publication-grade MCMC diagnostic plots (trace plots and posterior density plots)
#' for the spatial boundary changepoint and slope parameters, enabling rigorous Bayesian verification.
#' 
#' @param x A `spatialbp_fit` object.
#' @param warmup Integer, number of iterations to discard as warmup (default is 1000).
#' @param ... Unused.
#' @importFrom ggplot2 ggplot aes geom_line geom_density labs theme_minimal
#' @importFrom gridExtra grid.arrange
#' @export
plot_diagnostics <- function(x, warmup = 1000, ...) {
  UseMethod("plot_diagnostics")
}

#' @export
plot_diagnostics.spatialbp_fit <- function(x, warmup = 1000, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting diagnostics.")
  }
  
  fit_trace <- x$fit
  c_trace <- fit_trace$c
  ll_trace <- fit_trace$log_likelihood
  
  total_iters <- length(c_trace)
  burn_samples <- min(as.integer(warmup), total_iters - 10)
  
  # Filter out burn-in
  keep_indices <- (burn_samples + 1):total_iters
  c_trace_keep <- c_trace[keep_indices]
  ll_trace_keep <- ll_trace[keep_indices]
  
  # Create data frame for plotting
  df_trace <- data.frame(
    Iteration = keep_indices,
    c = c_trace_keep,
    LogLikelihood = ll_trace_keep
  )
  
  # 1. Trace plot of c
  p1 <- ggplot2::ggplot(df_trace, ggplot2::aes(x = Iteration, y = c)) +
    ggplot2::geom_line(color = "royalblue", alpha = 0.7) +
    ggplot2::labs(title = "Trace Plot of Boundary Changepoint (c)",
                  subtitle = paste("Warmup of", burn_samples, "samples discarded"),
                  x = "Iteration", y = "c") +
    ggplot2::theme_minimal()
  
  # 2. Posterior density of c
  p2 <- ggplot2::ggplot(df_trace, ggplot2::aes(x = c)) +
    ggplot2::geom_density(fill = "royalblue", alpha = 0.4, color = "royalblue") +
    ggplot2::labs(title = "Posterior Density of Changepoint (c)",
                  subtitle = paste("Warmup of", burn_samples, "samples discarded"),
                  x = "c", y = "Density") +
    ggplot2::theme_minimal()
  
  # 3. Log-likelihood trace
  p3 <- ggplot2::ggplot(df_trace, ggplot2::aes(x = Iteration, y = LogLikelihood)) +
    ggplot2::geom_line(color = "darkorange", alpha = 0.7) +
    ggplot2::labs(title = "Trace Plot of Log-Likelihood",
                  x = "Iteration", y = "Log-Likelihood") +
    ggplot2::theme_minimal()
  
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
  } else {
    message("gridExtra not installed. Returning list of individual plots.")
    return(list(trace_c = p1, density_c = p2, trace_ll = p3))
  }
}
