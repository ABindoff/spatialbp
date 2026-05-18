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
}
