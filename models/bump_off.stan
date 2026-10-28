// Constant-only model for a SINGLE trait -- the "bump off" alternative
// to bump_on.stan.
data {
  int<lower=1> N;
  int<lower=1> S;
  vector[N] y;
  array[N] int<lower=1, upper=S> site;
}

parameters {
  vector[S] delta;
  real mu_delta;
  real<lower=0> sd_delta;
  real<lower=0> sigma;
}

model {
  mu_delta ~ normal(0, 10);
  sd_delta ~ normal(0, 10);
  sigma    ~ normal(0, 10);

  delta ~ normal(mu_delta, sd_delta);
  y ~ normal(delta[site], sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = normal_lpdf(y[n] | delta[site[n]], sigma);
}
