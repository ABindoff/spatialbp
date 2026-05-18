# Spatial Affinity Conceptual Framework

## The Scenario
Imagine you send a satellite over a planet and imaging confirms three life forms which are identified by their colour: Red, Green, and Blue. Looking at the images it appears like Red and Green have an affinity, but Green and Blue do not like to co-habit an area of the planet. How would you test these hypotheses?

## Phase 1: Discrete Organisms (Point Processes)
If life forms are discrete entities, the framework relies on Multitype Spatial Point Patterns.

### Distance-Based Summary Statistics (Exploratory)
- **Affinity (Red & Green):** We expect the Cross L-function ($L_{RG}(r)$) to show a positive deviation from chance.
- **Repulsion (Green & Blue):** We expect the Cross L-function ($L_{GB}(r)$) to show a negative deviation (a "desert" effect).
- We test this via Monte Carlo simulations to create a 95% confidence envelope around the null hypothesis of independent spatial processes.

### Point Process Models (Inferential)
- Fit Multitype Gibbs Spatial Point Process Models to estimate interaction coefficients ($\beta_{RG} > 0$ for affinity, $\beta_{GB} < 0$ for repulsion).
- To avoid the "Hidden Variable" trap (where organisms share environmental preferences rather than actually interacting), we can use Joint Species Distribution Models (JSDMs) to model abundance against covariates while parsing out residual spatial correlations.

## Phase 2: Continuous Microorganism Colonies (Geostatistics)
If the colors represent continuous intensities of microorganism colonies rather than discrete points, the mathematical landscape shifts.

### Cross-Covariance Functions
- We replace Ripley's K with Cross-Variograms or Cross-Covariance Functions ($C(h)$).
- Affinity yields a strongly positive $C_{RG}(h)$ at short distances.
- Repulsion yields a strongly negative $C_{GB}(h)$ at short distances.

### Multivariate Gaussian Processes (MVGP)
- We model the intensities jointly using a Linear Model of Coregionalization (LMC).
- This yields full posterior distributions for the spatial range, cross-correlation strength, and the underlying fields.
- **The Bottleneck:** $\mathcal{O}(N^3)$ computational scaling.
- **The Solution:** Scalable Bayesian approximations like SPDE via INLA or Nearest Neighbor Gaussian Processes (NNGP).

## Phase 3: The `smoothbp()` Constraint
What if the only tool available on the satellite is `smoothbp()`?

### The Walls
1. **Dimensionality:** `smoothbp()` is 1D. Spatial fields are 2D. We must take 1D transect slices, destroying 2D autocorrelation.
2. **Breakpoints vs. Covariance:** `smoothbp()` models structural edge shifts, not continuous covariance. 

### The "Hacker" Workaround
- Draw hundreds of 1D transects.
- Run `smoothbp()` on Red, Green, and Blue independently per transect.
- **Test for Breakpoint Colocalization:**
  - *Affinity:* Compare posterior distributions of breakpoints. If `Breakpoint_Red == Breakpoint_Green` with high probability, they share borders.
  - *Repulsion:* A positive structural shift in Green predicts a negative structural shift in Blue.
