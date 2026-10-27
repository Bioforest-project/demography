// Bump model:  annual contribution to change = constant rate + gaussian bump
// inst_y = delta + alpha * exp(-0.5*((trec - tau)/gamma)^2) + Normal(0,sigma)

data {
  int<lower=1> N;                                // number of observations
  int<lower=1> S;                                // number of sites
  int<lower=1> P;                                // number of plots
  int<lower=1> K;                                // number of traits
  vector[N] y;                                   // response variable
  vector[N] trec;                                // recovery time
  vector[N] t0;                                  // recovery time of initial census
  array[N] int<lower=1, upper=S> site;           // site index for each observation
  array[N] int<lower=1, upper=P> plot;           // plot index for each observation
  array[N] int<lower=1, upper=K> trait;          // trait index for each observation
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
  matrix[S, K] delta;                            // baseline level per site x trait
  matrix[P, K] alpha;                            // bump amplitude per plot x trait
  vector<lower=tau_range[1], upper=tau_range[2]>[S] tau; // bump location (peak time) per site; cannot be outside of the measurement period
  vector<lower=gamma_range[1], upper=gamma_range[2]>[S] gamma;    // bump width (>0) per site; 95% CI cannot be smaller than 2 yrs (non detectable)
  vector<lower=0>[S] sigma;                      // residual sd per site (trait values are standardised)
}

model {
  vector[S * K] delta_flat = to_vector(delta);
  vector[P * K] alpha_flat = to_vector(alpha);

  // priors (weakly informative)
  delta_flat ~ normal(0, 1);
  alpha_flat ~ normal(0, 5);
  tau   ~ normal(10, 5);
  gamma ~ normal(3, 2);
  sigma ~ normal(0, 1);

  // likelihood
  vector[N] delta_n = delta_flat[idx_st];
  vector[N] alpha_n = alpha_flat[idx_pt];
  vector[N] tau_n   = tau[site];
  vector[N] gamma_n = gamma[site];
  vector[N] sigma_n = sigma[site];

  vector[N] mu = delta_n + alpha_n .* exp(-0.5 .* ((trec-tau_n)./gamma_n)^2) ;

  y ~ normal(mu, sigma_n);
}
