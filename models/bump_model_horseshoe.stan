// Bump-or-constant model with a regularized horseshoe prior (Piironen & Vehtari,
// 2017). Selection is at the trait level: each trait k gets a local shrinkage 
// scale lambda[k]; combined with a global scale tau0, this produces alpha[,k] 
// that either collapses toward 0 (trait doesn't need a bump) or is left free.

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> P;
  int<lower=1> K;
  vector[N] y;
  vector[N] trec;
  array[N] int<lower=1, upper=S> site;
  array[N] int<lower=1, upper=P> plot;
  array[N] int<lower=1, upper=K> trait;

  real tau_fixed;              // anchor value for tau when trait is "off"
  real<lower=0> sd_tau_off;    // sd around tau_fixed when "off"
  real<lower=0> gamma_fixed;   // anchor value for gamma when "off"
  real<lower=0> sd_gamma_off;  // sd around gamma_fixed when "off"

  // horseshoe hyperparameters (Piironen & Vehtari recommendations)
  real<lower=0> tau0_scale;    // global scale for tau0's half-Cauchy prior;
                                // smaller = more aggressive shrinkage.
  real<lower=0> slab_scale;    // c: slab scale (regularizes large alpha)
  real<lower=0> slab_df;       // nu: df for the slab-scale prior
}
transformed data {
  real max_trec = max(trec);
  array[N] int idx_st;
  array[N] int idx_pt;
  for (n in 1:N) {
    idx_st[n] = (trait[n] - 1) * S + site[n];
    idx_pt[n] = (trait[n] - 1) * P + plot[n];
  }
}
parameters {
  matrix[S, K] delta;
  // non-centered horseshoe parameterization for alpha
  matrix[P, K] alpha_tilde;     // raw N(0,1) draws
  vector<lower=0>[K] lambda;    // per-trait local scale (half-Cauchy)
  real<lower=0> tau0;           // global scale (half-Cauchy), shared across traits
  real<lower=0> c2;             // slab scale^2 (inverse-gamma)
  matrix<lower=0, upper=max_trec>[S, K] tau;
  vector<lower=0, upper=max_trec>[K] mu_tau;
  real<lower=0> sd_tau;
  matrix<lower=0>[S, K] gamma;
  vector<lower=0>[K] mu_gamma;
  real<lower=0> sd_gamma;
  vector<lower=0>[K] sigma;
}

transformed parameters {
  // regularized local scale: keeps alpha bounded even if tau0*lambda[k] is huge
  vector[K] lambda_tilde;
  for (k in 1:K)
    lambda_tilde[k] = sqrt(c2 * square(lambda[k]) / (c2 + square(tau0) * square(lambda[k])));
  // per-trait continuous "activation weight" in (0,1): near 0 -> trait off,
  // near 1 -> trait on. Replaces the discrete theta[k] indicator.
  vector[K] w;
  for (k in 1:K) {
    real s = tau0 * lambda_tilde[k];
    w[k] = s / (1 + s);
  }
  // alpha, non-centered: this parameterization (not alpha ~ normal(0, tau0*lambda))
  // is what makes horseshoe models sample efficiently
  matrix[P, K] alpha;
  for (k in 1:K)
    alpha[, k] = tau0 * lambda_tilde[k] * alpha_tilde[, k];
}
model {
  // --- horseshoe priors ---
  to_vector(alpha_tilde) ~ std_normal();
  lambda ~ cauchy(0, 1);
  tau0   ~ cauchy(0, tau0_scale);
  c2     ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df * square(slab_scale));
  // --- other priors ---
  to_vector(delta) ~ normal(0, 10);
  mu_tau   ~ normal(0, 50);
  sd_tau   ~ normal(0, 20);
  mu_gamma ~ normal(0, 20);
  sd_gamma ~ normal(0, 10);
  sigma    ~ normal(0, 10);
  
  // tau/gamma: continuous blend between "off anchor" and "on hierarchical",
  // driven by w[k] -- one smooth Normal density, no mixture.
  for (k in 1:K) {
    real tau_mean_k   = tau_fixed   + w[k] * (mu_tau[k]   - tau_fixed);
    real tau_sd_k      = sd_tau_off  + w[k] * (sd_tau       - sd_tau_off);
    real gamma_mean_k = gamma_fixed + w[k] * (mu_gamma[k] - gamma_fixed);
    real gamma_sd_k    = sd_gamma_off + w[k] * (sd_gamma     - sd_gamma_off);
    tau[, k]   ~ normal(tau_mean_k, tau_sd_k);
    gamma[, k] ~ normal(gamma_mean_k, gamma_sd_k);
  }
  // --- likelihood, vectorized ---
  {
    vector[S * K] delta_flat = to_vector(delta);
    vector[S * K] tau_flat   = to_vector(tau);
    vector[S * K] gamma_flat = to_vector(gamma);
    vector[P * K] alpha_flat = to_vector(alpha);

    vector[N] delta_n = delta_flat[idx_st];
    vector[N] tau_n   = tau_flat[idx_st];
    vector[N] gamma_n = gamma_flat[idx_st];
    vector[N] alpha_n = alpha_flat[idx_pt];
    vector[N] sigma_n = sigma[trait];

    vector[N] z  = (trec - tau_n) ./ gamma_n;
    vector[N] mu = delta_n + alpha_n .* exp(-0.5 * square(z));

    y ~ normal(mu, sigma_n);
  }
}
