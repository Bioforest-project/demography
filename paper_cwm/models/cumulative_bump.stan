// Bump model for cumulative sums, with residual autocorrelation modeled
// explicitly. Since cum_y is a running integral of (roughly) i.i.d. annual
// noise, the residual process within a plot is a Wiener process (Brownian
// motion) started at t0:
//   Var(resid_i)            = sigma^2 * (trec_i - t0)
//   Cov(resid_i, resid_j)   = sigma^2 * (min(trec_i, trec_j) - t0)
// for observations i, j from the same plot. 

functions {
  vector normal_cdf_vec(vector x, vector mu, vector sigma) {
    return 0.5 * (1 + erf((x - mu) ./ (sigma * sqrt(2.0))));
  }
}

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> P;
  int<lower=1> K;
  vector[N] cum_y;
  vector[N] trec;
  vector[N] t0;
  array[N] int<lower=1, upper=S> site;
  array[N] int<lower=1, upper=P> plot;
  array[N] int<lower=1, upper=K> trait;

  int<lower=1> B;                     // number of (plot, trait) blocks
  array[B] int<lower=1> block_start;  // 1-based start row of each block
  array[B] int<lower=1> block_size;   // number of rows in each block
  
  vector[2] tau_range;
  vector[2] gamma_range;
}

transformed data {
  array[N] int idx_st;
  array[N] int idx_pt;
  for (n in 1:N) {
    idx_st[n] = (trait[n] - 1) * S + site[n];
    idx_pt[n] = (trait[n] - 1) * P + plot[n];
  }
}

parameters {
  matrix[S, K] delta;
  matrix[P, K] alpha;
  vector<lower=tau_range[1], upper=tau_range[2]>[S] tau; // bump location (peak time) per site; cannot be outside of the measurement period
  vector<lower=gamma_range[1], upper=gamma_range[2]>[S] gamma;    // bump width (>0) per site; 95% CI cannot be smaller than 2 yrs (non detectable)
  matrix<lower=0>[P, K] sigma;         // residual sd per plot x trait
}

model {
  vector[S * K] delta_flat = to_vector(delta);
  vector[P * K] alpha_flat = to_vector(alpha);
  vector[P * K] sigma_flat = to_vector(sigma);

  delta_flat ~ normal(0, 5);
  alpha_flat ~ normal(0, 10);
  tau   ~ normal(10, 20);
  gamma ~ normal(0, 10);
  sigma_flat ~ normal(0, 5);

  // full-data vectorized mean function (still needed even with correlated
  // residuals -- only the covariance structure changes)
  vector[N] delta_n = delta_flat[idx_st];
  vector[N] alpha_n = alpha_flat[idx_pt];
  vector[N] tau_n   = tau[site];
  vector[N] gamma_n = gamma[site];
  vector[N] cdf_diff = normal_cdf_vec(trec, tau_n, gamma_n) - normal_cdf_vec(t0, tau_n, gamma_n);
  vector[N] mu = delta_n .* (trec - t0) + alpha_n .* gamma_n * sqrt(2 * pi()) .* cdf_diff;

  // --- likelihood, one multivariate normal per (plot, trait) block ---
  for (b in 1:B) {
    int start = block_start[b];
    int n_b   = block_size[b];
    int end   = start + n_b - 1;

    int p_b = plot[start];    // plot of this block (constant within block)
    int k_b = trait[start];   // trait of this block (constant within block)
    real sigma_b = sigma[p_b, k_b];

    vector[n_b] trec_b = trec[start:end];
    real t0_b = t0[start];    // assumed constant within block

    matrix[n_b, n_b] Sigma;
    for (i in 1:n_b) {
      for (j in 1:n_b) {
        Sigma[i, j] = square(sigma_b) * (fmin(trec_b[i], trec_b[j]) - t0_b);
      }
      Sigma[i, i] += 1e-8;   // small jitter for numerical stability
    }
    matrix[n_b, n_b] L_Sigma = cholesky_decompose(Sigma);

    cum_y[start:end] ~ multi_normal_cholesky(mu[start:end], L_Sigma);
  }
}