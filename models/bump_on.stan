// Bump model for a SINGLE trait (traits are independent in the likelihood,
// so fit one trait at a time and compare against bump_off.stan).
data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> P;
  vector[N] y;
  vector[N] trec;
  array[N] int<lower=1, upper=S> site;
  array[N] int<lower=1, upper=P> plot;
}

transformed data {
  real max_trec = max(trec);
}

parameters {
  vector[S] delta;
  real mu_delta;
  real<lower=0> sd_delta;

  vector[P] alpha;

  vector<lower=0, upper=max_trec>[S] tau;
  real<lower=0, upper=max_trec> mu_tau;
  real<lower=0> sd_tau;

  vector<lower=0>[S] gamma;
  real<lower=0> mu_gamma;
  real<lower=0> sd_gamma;

  real<lower=0> sigma;
}

model {
  mu_delta ~ normal(0, 10);
  sd_delta ~ normal(0, 10);
  mu_tau   ~ normal(0, 50);
  sd_tau   ~ normal(0, 20);
  mu_gamma ~ normal(0, 20);
  sd_gamma ~ normal(0, 10);
  sigma    ~ normal(0, 10);

  delta ~ normal(mu_delta, sd_delta);
  tau   ~ normal(mu_tau, sd_tau);
  gamma ~ normal(mu_gamma, sd_gamma);

  vector[N] z  = (trec - tau[site]) ./ gamma[site];
  vector[N] mu = delta[site] + alpha[plot] .* exp(-0.5 * square(z));
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  {
    vector[N] z  = (trec - tau[site]) ./ gamma[site];
    vector[N] mu = delta[site] + alpha[plot] .* exp(-0.5 * square(z));
    for (n in 1:N)
      log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
  }
}
