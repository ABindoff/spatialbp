use extendr_api::prelude::*;
use rand::distributions::{Distribution, Uniform};
use rand_distr::Normal;

/// Compute the Sum of Squared Errors (SSE) for the 2D Spatial Piecewise Model
/// @export
#[extendr]
fn compute_spatial_sse(
    observed: &[f64],
    baseline: &[f64],
    latent: &[f64],
    beta: f64,
    c: f64,
) -> f64 {
    let n = observed.len();
    // Prior check: c must lie within [min(latent), max(latent)]
    let mut min_z = f64::INFINITY;
    let mut max_z = f64::NEG_INFINITY;
    for &z in latent {
        if z < min_z { min_z = z; }
        if z > max_z { max_z = z; }
    }
    if c < min_z || c > max_z {
        return f64::INFINITY;
    }
    if baseline.len() != n || latent.len() != n {
        return f64::INFINITY;
    }

    let mut sse = 0.0;
    for i in 0..n {
        let hinge = if latent[i] > c { latent[i] - c } else { 0.0 };
        let y_hat = baseline[i] + beta * hinge;
        let error = observed[i] - y_hat;
        sse += error * error;
    }
    sse
}

/// Run a Random Walk Metropolis-Hastings MCMC to fit the spatial breakpoint
/// 
/// @param observed The observed pixel intensities (multiplexed histology data)
/// @param baseline The baseline spatial surface Z(x,y)
/// @param latent The latent propensity field
/// @param iterations Number of MCMC iterations
/// @param init_beta Starting value for beta
/// @param init_c Starting value for c
/// @param sd_beta Proposal standard deviation for beta
/// @param sd_c Proposal standard deviation for c
/// @param lambda L1 (Laplace) regularization parameter for beta
/// @export
#[extendr]
fn run_spatial_mcmc(
    observed: &[f64],
    baseline: &[f64],
    latent: &[f64],
    iterations: i32,
    init_beta: f64,
    init_c: f64,
    sd_beta: f64,
    sd_c: f64,
    lambda: f64,
) -> List {
    let mut rng = rand::thread_rng();
    let prop_beta_dist = Normal::new(0.0, sd_beta).unwrap();
    let prop_c_dist = Normal::new(0.0, sd_c).unwrap();
    let unif_dist = Uniform::new(0.0, 1.0);

    let mut current_beta = init_beta;
    let mut current_c = init_c;
    let current_sse = compute_spatial_sse(observed, baseline, latent, current_beta, current_c);
    
    // Log-Likelihood with Laplace Prior (L1 Regularization)
    let mut current_ll = -0.5 * current_sse - lambda * current_beta.abs();

    let iters = iterations as usize;
    let mut trace_beta = Vec::with_capacity(iters);
    let mut trace_c = Vec::with_capacity(iters);
    let mut trace_ll = Vec::with_capacity(iters);
    let mut acceptance_count = 0;

    for _ in 0..iters {
        // 1. Propose new state
        let prop_beta = current_beta + prop_beta_dist.sample(&mut rng);
        let prop_c = current_c + prop_c_dist.sample(&mut rng);

        // 2. Evaluate proposed likelihood
        let prop_sse = compute_spatial_sse(observed, baseline, latent, prop_beta, prop_c);
        let prop_ll = -0.5 * prop_sse - lambda * prop_beta.abs();

        // 3. Acceptance ratio
        let log_ratio = prop_ll - current_ll;
        
        let accept = if log_ratio > 0.0 {
            true
        } else {
            log_ratio.exp() > unif_dist.sample(&mut rng)
        };

        if accept {
            current_beta = prop_beta;
            current_c = prop_c;
            current_ll = prop_ll;
            acceptance_count += 1;
        }

        // 4. Record trace
        trace_beta.push(current_beta);
        trace_c.push(current_c);
        trace_ll.push(current_ll);
    }

    list!(
        beta = trace_beta,
        c = trace_c,
        log_likelihood = trace_ll,
        acceptance_rate = (acceptance_count as f64) / (iters as f64)
    )
}

fn compute_hinge_sse(
    observed: &[f64],
    baseline: &[f64],
    latent: &[f64],
    beta: f64,
    c: f64,
    dir: f64,
) -> f64 {
    let n = observed.len();
    // Prior check: c must lie within [min(latent), max(latent)]
    let mut min_z = f64::INFINITY;
    let mut max_z = f64::NEG_INFINITY;
    for &z in latent {
        if z < min_z { min_z = z; }
        if z > max_z { max_z = z; }
    }
    if c < min_z || c > max_z {
        return f64::INFINITY;
    }
    if baseline.len() != n || latent.len() != n {
        return f64::INFINITY;
    }

    let mut sse = 0.0;
    for i in 0..n {
        let val = dir * (latent[i] - c);
        let hinge = if val > 0.0 { val } else { 0.0 };
        let y_hat = baseline[i] + beta * hinge;
        let error = observed[i] - y_hat;
        sse += error * error;
    }
    sse
}

/// Fit a Joint 2D Spatial Piecewise Regression Model
/// 
/// Supports both "affinity" (shared boundary, same direction) and
/// "repulsion" (shared boundary, opposite directions) hypotheses.
/// 
/// @export
#[extendr]
fn run_joint_mcmc(
    y1: &[f64],
    y2: &[f64],
    baseline1: &[f64],
    baseline2: &[f64],
    latent: &[f64],
    mode: &str,
    iterations: i32,
    init_beta1: f64,
    init_beta2: f64,
    init_c: f64,
    sd_beta1: f64,
    sd_beta2: f64,
    sd_c: f64,
    lambda1: f64,
    lambda2: f64,
) -> List {
    let mut rng = rand::thread_rng();
    let prop_beta1_dist = Normal::new(0.0, sd_beta1).unwrap();
    let prop_beta2_dist = Normal::new(0.0, sd_beta2).unwrap();
    let prop_c_dist = Normal::new(0.0, sd_c).unwrap();
    let unif_dist = Uniform::new(0.0, 1.0);

    let mut current_beta1 = init_beta1;
    let mut current_beta2 = init_beta2;
    let mut current_c = init_c;

    let dir1 = 1.0;
    let dir2 = if mode == "repulsion" { -1.0 } else { 1.0 };

    let sse1 = compute_hinge_sse(y1, baseline1, latent, current_beta1, current_c, dir1);
    let sse2 = compute_hinge_sse(y2, baseline2, latent, current_beta2, current_c, dir2);
    
    // Log-Likelihood with Joint Laplace Priors (L1 Regularization)
    let mut current_ll = -0.5 * (sse1 + sse2) - lambda1 * current_beta1.abs() - lambda2 * current_beta2.abs();

    let iters = iterations as usize;
    let mut trace_beta1 = Vec::with_capacity(iters);
    let mut trace_beta2 = Vec::with_capacity(iters);
    let mut trace_c = Vec::with_capacity(iters);
    let mut trace_ll = Vec::with_capacity(iters);
    let mut acceptance_count = 0;

    for _ in 0..iters {
        // 1. Propose new state
        let prop_beta1 = current_beta1 + prop_beta1_dist.sample(&mut rng);
        let prop_beta2 = current_beta2 + prop_beta2_dist.sample(&mut rng);
        let prop_c = current_c + prop_c_dist.sample(&mut rng);

        // 2. Evaluate proposed likelihood
        let prop_sse1 = compute_hinge_sse(y1, baseline1, latent, prop_beta1, prop_c, dir1);
        let prop_sse2 = compute_hinge_sse(y2, baseline2, latent, prop_beta2, prop_c, dir2);
        
        let prop_ll = -0.5 * (prop_sse1 + prop_sse2) - lambda1 * prop_beta1.abs() - lambda2 * prop_beta2.abs();

        // 3. Acceptance ratio
        let log_ratio = prop_ll - current_ll;
        let accept = if log_ratio > 0.0 {
            true
        } else {
            log_ratio.exp() > unif_dist.sample(&mut rng)
        };

        if accept {
            current_beta1 = prop_beta1;
            current_beta2 = prop_beta2;
            current_c = prop_c;
            current_ll = prop_ll;
            acceptance_count += 1;
        }

        // 4. Record trace
        trace_beta1.push(current_beta1);
        trace_beta2.push(current_beta2);
        trace_c.push(current_c);
        trace_ll.push(current_ll);
    }

    list!(
        beta1 = trace_beta1,
        beta2 = trace_beta2,
        c = trace_c,
        log_likelihood = trace_ll,
        acceptance_rate = (acceptance_count as f64) / (iters as f64)
    )
}

extendr_module! {
    mod spatialbp;
    fn compute_spatial_sse;
    fn run_spatial_mcmc;
    fn run_joint_mcmc;
    fn run_formula_mcmc;
}
/// Compute SSE for the extended Formula model
fn compute_formula_sse(
    y: &[f64],
    offset: &[f64],
    baseline_matrix: &[f64],
    latent: &[f64],
    beta_covariates: &[f64],
    c_covariates: &[f64],
    baseline_w: &[f64],
    beta_w: &[f64],
    c_w: &[f64],
    smooth: bool,
    rho: f64,
    dir: f64,
    k: usize,
) -> f64 {
    let n = latent.len();
    let p = baseline_w.len() / k;
    let b = beta_w.len() / k;
    let c_dim = c_w.len();
    
    let mut sse = 0.0;
    
    for i in 0..n {
        // 1. Calculate localized c (Shared across all K dimensions)
        let mut c_i = 0.0;
        for j in 0..c_dim {
            c_i += c_covariates[i * c_dim + j] * c_w[j];
        }
        
        // Evaluate Hinge (Shared boundary)
        let val = dir * (latent[i] - c_i);
        let hinge = if smooth {
            if rho <= 0.0001 {
                if val > 0.0 { val } else { 0.0 }
            } else {
                let scaled = val / rho;
                if scaled > 50.0 { val } else { rho * (1.0 + scaled.exp()).ln() }
            }
        } else {
            if val > 0.0 { val } else { 0.0 }
        };
        
        // Iterate over each response dimension k
        for dim in 0..k {
            let mut beta_i = 0.0;
            for j in 0..b {
                beta_i += beta_covariates[i * b + j] * beta_w[dim * b + j];
            }
            
            let mut base_i = offset[i + dim * n];
            for j in 0..p {
                base_i += baseline_matrix[i * p + j] * baseline_w[dim * p + j];
            }
            
            let y_hat = base_i + beta_i * hinge;
            let error = y[i + dim * n] - y_hat;
            sse += error * error;
        }
    }
    
    sse
}

#[extendr]
pub fn run_formula_mcmc(
    y: &[f64],
    offset: &[f64],
    baseline_matrix: &[f64],
    latent: &[f64],
    beta_covariates: &[f64],
    c_covariates: &[f64],
    num_baseline: i32,
    num_beta: i32,
    num_c: i32,
    num_responses: i32,
    smooth: bool,
    iterations: i32,
    lambda: f64,
    sd_beta: f64,
    sd_c: f64,
) -> List {
    let mut rng = rand::thread_rng();
    let unif_dist = Uniform::new(0.0, 1.0);
    let prop_beta_dist = Normal::new(0.0, sd_beta).unwrap();
    let prop_c_dist = Normal::new(0.0, sd_c).unwrap();
    let prop_rho_dist = Normal::new(0.0, 0.05).unwrap();
    // Default tiny step size for baseline to avoid destabilizing the chain
    let prop_baseline_dist = Normal::new(0.0, 0.05).unwrap();

    let k = num_responses as usize;
    let p = num_baseline as usize;
    let b = num_beta as usize;
    let c_dim = num_c as usize;

    let mut current_baseline = vec![0.0; p * k];
    let mut current_beta = vec![0.0; b * k];
    
    let n = latent.len();
    let n_f64 = n as f64;
    for dim in 0..k {
        let mut sum_y = 0.0;
        let mut sum_z = 0.0;
        for i in 0..n {
            sum_y += y[i + dim * n];
            sum_z += latent[i];
        }
        let mean_y = sum_y / n_f64;
        let mean_z = sum_z / n_f64;
        
        let mut cov_yz = 0.0;
        let mut var_z = 0.0;
        for i in 0..n {
            let dz = latent[i] - mean_z;
            cov_yz += (y[i + dim * n] - mean_y) * dz;
            var_z += dz * dz;
        }
        
        let mut init_beta = if var_z > 1e-6 { cov_yz / var_z } else { 0.0 };
        // We expect positive affinity in most user cases, but OLS gives the best start.
        if b > 0 {
            current_beta[dim * b] = init_beta;
        }
        if p > 0 {
            // Initialize intercept assuming hinge is roughly zero-centered
            current_baseline[dim * p] = mean_y;
        }
    }
    
    let mut current_c = vec![0.0; c_dim];
    
    let mut min_z = f64::INFINITY;
    let mut max_z = f64::NEG_INFINITY;
    
    if c_dim > 0 && latent.len() > 0 {
        let mut sum_latent = 0.0;
        for &z in latent {
            if z < min_z { min_z = z; }
            if z > max_z { max_z = z; }
            sum_latent += z;
        }
        current_c[0] = sum_latent / latent.len() as f64;
        
        // Ensure initial c is within absolute limits
        if current_c[0] < min_z { current_c[0] = min_z; }
        if current_c[0] > max_z { current_c[0] = max_z; }
    }

    let mut current_rho = if smooth { 1.0 } else { 0.0001 };

    let mut current_sse = compute_formula_sse(
        y, offset, baseline_matrix, latent, beta_covariates, c_covariates,
        &current_baseline, &current_beta, &current_c, smooth, current_rho, 1.0, k
    );
    let mut current_ll = -0.5 * current_sse;

    let mut acceptance_count = 0;
    
    let mut trace_beta0 = Vec::with_capacity(iterations as usize);
    let mut trace_c0 = Vec::with_capacity(iterations as usize);
    let mut trace_rho = Vec::with_capacity(iterations as usize);
    let mut trace_ll = Vec::with_capacity(iterations as usize);

    for _ in 0..iterations {
        // Block 1: Baseline
        if p > 0 {
            let mut prop_baseline = current_baseline.clone();
            for i in 0..(p * k) { prop_baseline[i] += prop_baseline_dist.sample(&mut rng); }
            let prop_sse = compute_formula_sse(
                y, offset, baseline_matrix, latent, beta_covariates, c_covariates,
                &prop_baseline, &current_beta, &current_c, smooth, current_rho, 1.0, k
            );
            let mut penalty = 0.0;
            for dim in 0..k { for j in 0..b { penalty += lambda * current_beta[dim * b + j].abs(); } }
            let prop_ll = -0.5 * prop_sse - penalty;
            let log_ratio = prop_ll - current_ll;
            if log_ratio > 0.0 || log_ratio.exp() > unif_dist.sample(&mut rng) {
                current_baseline = prop_baseline;
                current_ll = prop_ll;
                acceptance_count += 1;
            }
        }

        // Block 2: Beta
        if b > 0 {
            let mut prop_beta = current_beta.clone();
            for i in 0..(b * k) { prop_beta[i] += prop_beta_dist.sample(&mut rng); }
            let prop_sse = compute_formula_sse(
                y, offset, baseline_matrix, latent, beta_covariates, c_covariates,
                &current_baseline, &prop_beta, &current_c, smooth, current_rho, 1.0, k
            );
            let mut penalty = 0.0;
            for dim in 0..k { for j in 0..b { penalty += lambda * prop_beta[dim * b + j].abs(); } }
            let prop_ll = -0.5 * prop_sse - penalty;
            let log_ratio = prop_ll - current_ll;
            if log_ratio > 0.0 || log_ratio.exp() > unif_dist.sample(&mut rng) {
                current_beta = prop_beta;
                current_ll = prop_ll;
                acceptance_count += 1;
            }
        }

        // Block 3: C
        if c_dim > 0 {
            let mut prop_c = current_c.clone();
            for i in 0..c_dim { prop_c[i] += prop_c_dist.sample(&mut rng); }
            let prop_ll = if prop_c[0] < min_z || prop_c[0] > max_z {
                -f64::INFINITY
            } else {
                let prop_sse = compute_formula_sse(
                    y, offset, baseline_matrix, latent, beta_covariates, c_covariates,
                    &current_baseline, &current_beta, &prop_c, smooth, current_rho, 1.0, k
                );
                let mut penalty = 0.0;
                for dim in 0..k { for j in 0..b { penalty += lambda * current_beta[dim * b + j].abs(); } }
                -0.5 * prop_sse - penalty
            };
            let log_ratio = prop_ll - current_ll;
            if log_ratio > 0.0 || log_ratio.exp() > unif_dist.sample(&mut rng) {
                current_c = prop_c;
                current_ll = prop_ll;
                acceptance_count += 1;
            }
        }

        // Block 4: Rho
        if smooth {
            let mut prop_rho = current_rho + prop_rho_dist.sample(&mut rng);
            if prop_rho < 0.0001 { prop_rho = 0.0001; }
            let prop_sse = compute_formula_sse(
                y, offset, baseline_matrix, latent, beta_covariates, c_covariates,
                &current_baseline, &current_beta, &current_c, smooth, prop_rho, 1.0, k
            );
            let mut penalty = 0.0;
            for dim in 0..k { for j in 0..b { penalty += lambda * current_beta[dim * b + j].abs(); } }
            let prop_ll = -0.5 * prop_sse - penalty;
            let log_ratio = prop_ll - current_ll;
            if log_ratio > 0.0 || log_ratio.exp() > unif_dist.sample(&mut rng) {
                current_rho = prop_rho;
                current_ll = prop_ll;
                acceptance_count += 1;
            }
        }

        if b > 0 { trace_beta0.push(current_beta[0]); }
        if c_dim > 0 { trace_c0.push(current_c[0]); }
        if smooth { trace_rho.push(current_rho); }
        trace_ll.push(current_ll);
    }

    list!(
        beta0 = trace_beta0,
        c0 = trace_c0,
        rho = trace_rho,
        log_likelihood = trace_ll,
        acceptance_rate = (acceptance_count as f64) / (iterations as f64)
    )
}
