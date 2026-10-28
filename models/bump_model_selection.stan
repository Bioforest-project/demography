// Bump-or-constant model, site-level delta/tau/gamma, plot-level alpha.
// Spike-and-slab selection shared (via theta[k]) across alpha, tau, and gamma:
// when the bump is "off" for trait k, alpha shrinks toward 0 AND tau/gamma
// are pulled toward fixed anchor values instead of drifting freely under
// their hierarchical priors (which are unidentified when alpha ~ 0).

functions {
  // Elementwise log_mix over vectors, using the numerically-stable identity
  // log_sum_exp(a, b) = b + log1p_exp(a - b). Equivalent to calling
  // log_mix(theta, lp1[i], lp2[i]) for each i, but vectorized.
  vector log_mix_vec(real theta, vector lp1, vector lp2) {
    vector[num_elements(lp1)] a = log(theta) + lp1;
    vector[num_elements(lp1)] b = log1m(theta) + lp2;
    return b + log1p_exp(a - b);
  }
}

data {
  int<lower=1> N;
  int<lower=1> S;                              // number of sites
  int<lower=1> P;                              // number of plots
  int<lower=1> K;                              // number of traits
  vector[N] y;
  vector[N] trec;
  array[N] int<lower=1, upper=S> site;
  array[N] int<lower=1, upper=P> plot;
  array[N] int<lower=1, upper=K> trait;

  real<lower=0> sd_spike;    // "on"  width for alpha
  real<lower=0> sd_slab;   // "off" width for alpha

  real tau_fixed;              // anchor value for tau when bump is "off"
  real<lower=0> sd_tau_off;    // tight sd around tau_fixed when "off"
  real<lower=0> gamma_fixed;   // anchor value for gamma when bump is "off"
  real<lower=0> sd_gamma_off;  // tight sd around gamma_fixed when "off"
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
  real mu_delta;
  real<lower=0> sd_delta;
  matrix[P, K] alpha;
  matrix<lower=0, upper=max_trec>[S, K] tau;
  real<lower=0, upper=max_trec> mu_tau;
  real<lower=0> sd_tau;
  matrix<lower=0>[S, K] gamma;
  real<lower=0> mu_gamma;
  real<lower=0> sd_gamma;
  vector<lower=0>[K] sigma;
  vector<lower=0, upper=1>[K] theta;   // per-trait probability bump is "on"
}
model {
  // --- priors ---
  mu_delta ~ normal(0, 10);
  sd_delta ~ normal(0, 10);
  mu_tau   ~ normal(10, 20);
  sd_tau   ~ normal(0, 10);
  mu_gamma ~ normal(0, 10);
  sd_gamma ~ normal(0, 5);
  sigma    ~ normal(0, 10);
  theta    ~ beta(1, 1);

  // hierarchical model for delta
  to_vector(delta) ~ normal(mu_delta, sd_delta);

  // --- spike-and-slab, shared theta[k], applied to alpha, tau, gamma ---
  // (constants dropped from the normal log-density: they're identical in
  // both mixture components so they cancel in the relative weighting)
  for (k in 1:K) {
    vector[P] lp1_a = -log(sd_spike)  - 0.5 * square(alpha[, k] / sd_spike);
    vector[P] lp2_a = -log(sd_slab) - 0.5 * square(alpha[, k] / sd_slab);
    target += sum(log_mix_vec(theta[k], lp1_a, lp2_a));

    vector[S] lp1_t = -log(sd_tau)     - 0.5 * square((tau[, k]   - mu_tau)  / sd_tau);
    vector[S] lp2_t = -log(sd_tau_off) - 0.5 * square((tau[, k]   - tau_fixed)  / sd_tau_off);
    target += sum(log_mix_vec(theta[k], lp1_t, lp2_t));

    vector[S] lp1_g = -log(sd_gamma)     - 0.5 * square((gamma[, k] - mu_gamma) / sd_gamma);
    vector[S] lp2_g = -log(sd_gamma_off) - 0.5 * square((gamma[, k] - gamma_fixed) / sd_gamma_off);
    target += sum(log_mix_vec(theta[k], lp1_g, lp2_g));
  }

  // --- likelihood ---
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

generated quantities {
  // posterior probability that plot/site p,k has the bump "on"
  matrix[P, K] prob_bump_on;
  for (k in 1:K) {
    for (p in 1:P) {
      real lp_on  = log(theta[k])   - log(sd_spike)  - 0.5 * square(alpha[p, k] / sd_spike);
      real lp_off = log1m(theta[k]) - log(sd_slab) - 0.5 * square(alpha[p, k] / sd_slab);
      prob_bump_on[p, k] = exp(lp_on - log_sum_exp(lp_on, lp_off));
    }
  }
}
