#' Multivariate Spatial Affinity Test
#' 
#' Tests the affinity and repulsion hypotheses between three protein stains
#' (Red, Green, Blue) using the spatialbp Rust MCMC engine.
#' 
#' @param red Numeric vector of Red protein intensities
#' @param green Numeric vector of Green protein intensities
#' @param blue Numeric vector of Blue protein intensities
#' @param latent_field Precomputed spatial radial basis values
#' @param iterations Number of MCMC iterations
#' @return A list containing Bayes Factors and posterior probabilities.
#' @export
test_spatialbp <- function(red, green, blue, latent_field, iterations = 10000) {
  
  # --- MODEL 0: The Independent Null Hypothesis ---
  # Assuming no affinity. Each protein has its own independent spatial boundary.
  # We run the Rust MCMC separately for each to get the baseline log-likelihood.
  message("Fitting Independent Null Model...")
  mcmc_r_indep <- run_spatial_mcmc(red, rep(0, length(red)), latent_field, iterations, 1, 0, 0.1, 0.1)
  mcmc_g_indep <- run_spatial_mcmc(green, rep(0, length(green)), latent_field, iterations, 1, 0, 0.1, 0.1)
  mcmc_b_indep <- run_spatial_mcmc(blue, rep(0, length(blue)), latent_field, iterations, 1, 0, 0.1, 0.1)
  
  ll_indep <- mean(mcmc_r_indep$log_likelihood) + 
              mean(mcmc_g_indep$log_likelihood) + 
              mean(mcmc_b_indep$log_likelihood)
  
  # --- MODEL 1: The Affinity Hypothesis (Red & Green Share a Boundary) ---
  # We hypothesize that Red and Green share the EXACT same spatial contour and direction.
  message("Fitting Joint Affinity Model (Red & Green)...")
  mcmc_affinity <- run_joint_mcmc(
    y1 = red, y2 = green,
    baseline1 = rep(0, length(red)), baseline2 = rep(0, length(green)),
    latent = latent_field,
    mode = "affinity",
    iterations = iterations,
    init_beta1 = 1, init_beta2 = 1, init_c = 0,
    sd_beta1 = 0.1, sd_beta2 = 0.1, sd_c = 0.1
  )
  # Joint log-likelihood for Red + Green + Independent Blue
  ll_affinity <- mean(mcmc_affinity$log_likelihood) + mean(mcmc_b_indep$log_likelihood)
  
  # --- MODEL 2: The Repulsion Hypothesis (Green & Blue Exclude) ---
  # We hypothesize that Green and Blue share a boundary but occupy opposite sides of it.
  message("Fitting Joint Repulsion Model (Green & Blue)...")
  mcmc_repulsion <- run_joint_mcmc(
    y1 = green, y2 = blue,
    baseline1 = rep(0, length(green)), baseline2 = rep(0, length(blue)),
    latent = latent_field,
    mode = "repulsion",
    iterations = iterations,
    init_beta1 = 1, init_beta2 = 1, init_c = 0,
    sd_beta1 = 0.1, sd_beta2 = 0.1, sd_c = 0.1
  )
  # Joint log-likelihood for Independent Red + Joint Green/Blue
  ll_repulsion <- mean(mcmc_r_indep$log_likelihood) + mean(mcmc_repulsion$log_likelihood)
  
  # Calculate Bayes Factors (comparing joint models vs independent null models)
  # In Bayesian model comparison, Bayes Factor = exp(mean_log_lik_joint - mean_log_lik_indep)
  # (Approximated here via harmonic mean of likelihoods or posterior mean of log-likelihoods as a proxy)
  bayes_factor_affinity <- exp(ll_affinity - ll_indep)
  bayes_factor_repulsion <- exp(ll_repulsion - ll_indep)
  
  list(
    ll_null = ll_indep,
    ll_affinity = ll_affinity,
    ll_repulsion = ll_repulsion,
    bayes_factor_affinity = bayes_factor_affinity,
    bayes_factor_repulsion = bayes_factor_repulsion,
    mcmc_affinity = mcmc_affinity,
    mcmc_repulsion = mcmc_repulsion
  )
}

#' Automatically Discover Pairwise Spatial Relationships
#' 
#' Performs pairwise joint MCMC modeling across all provided protein stains to
#' automatically discover which pairs exhibit Affinity, Repulsion, or Independence.
#' 
#' @param stains A named list or data.frame where each element is a numeric vector of intensities.
#' @param latent_field Precomputed spatial radial basis or coordinate values.
#' @param iterations Number of MCMC iterations (default = 5000).
#' @param bf_threshold Numeric, Bayes Factor threshold to classify a relationship as significant (default = 10).
#' @return A list containing a relationship data.frame and a pairwise Bayes Factor matrix.
#' @export
discover_spatial_relationships <- function(stains, latent_field, iterations = 5000, bf_threshold = 10) {
  stain_names <- names(stains)
  n_stains <- length(stains)
  
  if (n_stains < 2) {
    stop("Need at least two stains to discover relationships.")
  }
  
  # 1. Precompute Independent Fits to establish baseline log-likelihoods
  message("Precomputing independent baseline models...")
  indep_ll <- list()
  for (name in stain_names) {
    message("  Fitting independent model for: ", name)
    fit_indep <- run_spatial_mcmc(
      observed = stains[[name]],
      baseline = rep(0, length(stains[[name]])),
      latent = latent_field,
      iterations = iterations,
      init_beta = 1.0, init_c = 0.0,
      sd_beta = 0.1, sd_c = 0.1, lambda = 0.0
    )
    indep_ll[[name]] <- mean(fit_indep$log_likelihood)
  }
  
  # Initialize results structures
  relationships <- data.frame(
    Stain1 = character(),
    Stain2 = character(),
    BF_Affinity = numeric(),
    BF_Repulsion = numeric(),
    Classification = character(),
    stringsAsFactors = FALSE
  )
  
  bf_matrix_affinity <- matrix(1, nrow = n_stains, ncol = n_stains, dimnames = list(stain_names, stain_names))
  bf_matrix_repulsion <- matrix(1, nrow = n_stains, ncol = n_stains, dimnames = list(stain_names, stain_names))
  
  # 2. Pairwise MCMC runs
  message("\nStarting pairwise relationship discovery...")
  pairs <- combn(stain_names, 2, simplify = FALSE)
  
  for (pair in pairs) {
    s1 <- pair[1]
    s2 <- pair[2]
    message("\nAnalyzing pair: ", s1, " <-> ", s2)
    
    # Baseline joint null (independent sum)
    ll_null <- indep_ll[[s1]] + indep_ll[[s2]]
    
    # Fit Affinity
    message("  Fitting Joint Affinity...")
    fit_aff <- run_joint_mcmc(
      y1 = stains[[s1]], y2 = stains[[s2]],
      baseline1 = rep(0, length(stains[[s1]])), baseline2 = rep(0, length(stains[[s2]])),
      latent = latent_field,
      mode = "affinity",
      iterations = iterations,
      init_beta1 = 1, init_beta2 = 1, init_c = 0,
      sd_beta1 = 0.1, sd_beta2 = 0.1, sd_c = 0.1,
      lambda1 = 0.0, lambda2 = 0.0
    )
    ll_affinity <- mean(fit_aff$log_likelihood)
    bf_aff <- exp(ll_affinity - ll_null)
    
    # Fit Repulsion
    message("  Fitting Joint Repulsion...")
    fit_rep <- run_joint_mcmc(
      y1 = stains[[s1]], y2 = stains[[s2]],
      baseline1 = rep(0, length(stains[[s1]])), baseline2 = rep(0, length(stains[[s2]])),
      latent = latent_field,
      mode = "repulsion",
      iterations = iterations,
      init_beta1 = 1, init_beta2 = 1, init_c = 0,
      sd_beta1 = 0.1, sd_beta2 = 0.1, sd_c = 0.1,
      lambda1 = 0.0, lambda2 = 0.0
    )
    ll_repulsion <- mean(fit_rep$log_likelihood)
    bf_rep <- exp(ll_repulsion - ll_null)
    
    # Save to matrices
    idx1 <- which(stain_names == s1)
    idx2 <- which(stain_names == s2)
    bf_matrix_affinity[idx1, idx2] <- bf_aff
    bf_matrix_affinity[idx2, idx1] <- bf_aff
    bf_matrix_repulsion[idx1, idx2] <- bf_rep
    bf_matrix_repulsion[idx2, idx1] <- bf_rep
    
    # Classification logic based on Bayes Factors
    classification <- "Independent"
    if (bf_aff > bf_threshold && bf_aff > bf_rep) {
      classification <- "Affinity"
    } else if (bf_rep > bf_threshold && bf_rep > bf_aff) {
      classification <- "Repulsion"
    }
    
    relationships <- rbind(relationships, data.frame(
      Stain1 = s1,
      Stain2 = s2,
      BF_Affinity = bf_aff,
      BF_Repulsion = bf_rep,
      Classification = classification,
      stringsAsFactors = FALSE
    ))
  }
  
  return(list(
    relationships = relationships,
    bf_matrix_affinity = bf_matrix_affinity,
    bf_matrix_repulsion = bf_matrix_repulsion
  ))
}
