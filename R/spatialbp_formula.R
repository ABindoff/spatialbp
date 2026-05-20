#' Define a spatial hinge term for the formula interface
#' 
#' @param latent The latent spatial gradient variable
#' @param mod_beta One-sided formula specifying covariates that modify the slope
#' @param mod_c One-sided formula specifying covariates that modify the breakpoint
#' @param smooth Logical. If TRUE, estimates a smooth transition (Softplus rho)
#' @export
hinge <- function(latent, mod_beta = NULL, mod_c = NULL, smooth = FALSE) {
  list(
    latent = latent,
    mod_beta = mod_beta,
    mod_c = mod_c,
    smooth = smooth
  )
}

#' Fit a Spatial Piecewise Model via Formula
#' 
#' @param formula A formula object (e.g., y ~ x + offset(z) + hinge(latent_z))
#' @param data A data frame
#' @param mode Character, "independent", "affinity", "repulsion"
#' @param iterations Integer, number of MCMC iterations
#' @export
spatialbp <- function(formula, data, mode = "independent", iterations = 5000, warmup = floor(iterations / 2), sd_beta = 0.1, sd_c = NULL) {
  message("Parsing spatialbp formula...")
  
  # 1. Extract special terms
  trm <- terms(formula, specials = c("hinge", "offset"))
  hinge_idx <- attr(trm, "specials")$hinge
  offset_idx <- attr(trm, "specials")$offset
  
  if (is.null(hinge_idx)) {
    stop("Formula must contain a hinge() term. Example: y ~ 1 + hinge(latent_z)")
  }
  
  if (length(hinge_idx) > 1) {
    stop("spatialbp currently supports exactly one hinge() term per formula.")
  }
  
  # 2. Extract Response(s)
  # Handle joint models like cbind(Green, Red)
  response_name <- rownames(attr(trm, "factors"))[1]
  y <- eval(parse(text = response_name), envir = data)
  
  # Ensure y is a matrix for consistency (multivariate or univariate)
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  
  # 3. Evaluate the hinge term
  hinge_call <- attr(trm, "variables")[[hinge_idx + 1]]
  hinge_eval <- eval(hinge_call, envir = data, enclos = parent.frame())
  
  latent_field <- hinge_eval$latent
  is_smooth <- hinge_eval$smooth
  
  # Parse slope modifiers
  beta_covariates <- matrix(1, nrow = nrow(data), ncol = 1) # Intercept only
  if (!is.null(hinge_eval$mod_beta)) {
    mf_beta <- model.frame(hinge_eval$mod_beta, data = data)
    beta_covariates <- model.matrix(hinge_eval$mod_beta, mf_beta)
  }
  
  # Parse breakpoint modifiers
  c_covariates <- matrix(1, nrow = nrow(data), ncol = 1) # Intercept only
  if (!is.null(hinge_eval$mod_c)) {
    mf_c <- model.frame(hinge_eval$mod_c, data = data)
    c_covariates <- model.matrix(hinge_eval$mod_c, mf_c)
  }
  
  # 4. Extract Baseline Covariates (Everything except hinge)
  # We construct a new formula without the hinge term
  hinge_str <- paste(deparse(hinge_call), collapse = " ")
  baseline_formula <- update(formula, paste(". ~ . -", hinge_str))
  mf_baseline <- model.frame(baseline_formula, data, na.action = na.pass)
  
  # Extract standard model matrix for baseline
  baseline_matrix <- model.matrix(baseline_formula, mf_baseline)
  
  # Extract offset
  offset_val <- model.offset(mf_baseline)
  if (is.null(offset_val)) {
    offset_val <- matrix(0, nrow = nrow(data), ncol = ncol(y))
  } else {
    offset_val <- as.matrix(offset_val)
    if (ncol(offset_val) == 1 && ncol(y) > 1) {
      # Broadcast offset to all responses
      offset_val <- matrix(offset_val, nrow = nrow(data), ncol = ncol(y))
    }
  }
  
  message("Formula parsed successfully!")
  message(" - Responses: ", ncol(y))
  message(" - Baseline Covariates: ", ncol(baseline_matrix) - 1) # Excluding intercept
  message(" - Slope Modifiers: ", ncol(beta_covariates) - 1)
  message(" - Boundary Modifiers: ", ncol(c_covariates) - 1)
  message(" - Smooth Transition (rho): ", is_smooth)
  
  # ==========================================
  # ==========================================
  # Call Rust Engine (Supports K >= 1 Responses)
  # ==========================================
  message(sprintf("Calling Rust MCMC Backend for %d joint responses...", ncol(y)))
  
  # Scale Responses to ensure L1 Penalty (lambda) has mathematically consistent weight
  y_means <- apply(y, 2, mean)
  y_sds <- apply(y, 2, sd)
  y_sds[y_sds == 0] <- 1.0
  y_scaled <- y
  for (dim in seq_len(ncol(y))) {
    y_scaled[, dim] <- (y[, dim] - y_means[dim]) / y_sds[dim]
  }
  
  if (is.null(sd_c)) {
    sd_c <- sd(latent_field) * 0.05
    if (is.na(sd_c) || sd_c == 0) sd_c <- 0.1
  }

  fit <- run_formula_mcmc(
    y = as.numeric(y_scaled),
    offset = as.numeric(offset_val),
    baseline_matrix = as.numeric(baseline_matrix),
    latent = as.numeric(latent_field),
    beta_covariates = as.numeric(beta_covariates),
    c_covariates = as.numeric(c_covariates),
    num_baseline = as.integer(ncol(baseline_matrix)),
    num_beta = as.integer(ncol(beta_covariates)),
    num_c = as.integer(ncol(c_covariates)),
    num_responses = as.integer(ncol(y)),
    smooth = is_smooth,
    iterations = as.integer(iterations),
    lambda = 1.0, # Gentle L1 penalty
    sd_beta = as.numeric(sd_beta),
    sd_c = as.numeric(sd_c)
  )
  
  # Format for S3 Inheritance (smoothbp)
  formatted_fit <- list(
    traces = list(
      beta = fit$beta0 * y_sds[1], # Unscale primary slope for the first dimension
      cp = fit$c0
    ),
    log_likelihood = fit$log_likelihood,
    acceptance_rate = fit$acceptance_rate,
    warmup = as.integer(warmup),
    data = data,
    spatial_context = list(
      observed = y,
      latent = latent_field
    ),
    rho_trace = if (is_smooth) fit$rho else NULL
  )
  class(formatted_fit) <- c("spatialbp", "smoothbp")
  return(formatted_fit)
}

#' Confidence Intervals for spatialbp Model
#' 
#' @param object A fitted spatialbp model
#' @param parm Which parameter to compute confidence intervals for (currently supports "c")
#' @param level Confidence level (default 0.95)
#' @param ... Additional arguments
#' @export
confint.spatialbp <- function(object, parm = "c", level = 0.95, ...) {
  if (parm == "c") {
    # Extract the trace after warmup
    traces <- object$traces$cp
    if (length(traces) > object$warmup) {
      traces <- traces[(object$warmup + 1):length(traces)]
    }
    alpha <- 1 - level
    bounds <- quantile(traces, probs = c(alpha/2, 1 - alpha/2))
    names(bounds) <- c(paste0(alpha/2 * 100, "%"), paste0((1 - alpha/2) * 100, "%"))
    return(bounds)
  } else {
    stop("Only 'c' parameter is supported for spatialbp confint currently.")
  }
}

#' Compute Distance Space
#' 
#' Computes the Euclidean distance from every point (x, y) to a given center (cx, cy).
#' 
#' @param x Numeric vector of x coordinates
#' @param y Numeric vector of y coordinates
#' @param cx X coordinate of the center
#' @param cy Y coordinate of the center
#' @return Numeric vector of distances
#' @export
distance_space <- function(x, y, cx, cy) {
  sqrt((x - cx)^2 + (y - cy)^2)
}

#' Compute Density Space
#' 
#' Computes a simple Gaussian-kernel smoothed density space based on point intensities.
#' 
#' @param x Numeric vector of x coordinates
#' @param y Numeric vector of y coordinates
#' @param z Numeric vector of intensities or weights
#' @param bandwidth Numeric, the standard deviation of the Gaussian kernel
#' @return Numeric vector of smoothed densities
#' @export
density_space <- function(x, y, z, bandwidth = 5) {
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    d2 <- (x - x[i])^2 + (y - y[i])^2
    w <- exp(-d2 / (2 * bandwidth^2))
    out[i] <- sum(w * z)
  }
  out
}
