data {
  int<lower=0> n; // number of old-growth observations
  int<lower=0> s;
  vector[n] stocks;
  vector[n] influx;
  vector[n] outflux;
  array[n] int<lower=0, upper=s> site;   // site index
  array[2] real mu_theta_bounds;
  array[2] real theta_bounds;
  array[2] real omega_bounds;
}
parameters {
  real<lower=mu_theta_bounds[1], upper=mu_theta_bounds[2]> mu_theta; 
  real<lower=0> sigma_theta;
  vector<lower=theta_bounds[1], upper=theta_bounds[2]>[s] theta_s;
  real<lower=0> sigma_stocks;
  real<lower=0> sigma_influx;
  real<lower=0> sigma_outflux;
}
transformed parameters {
  vector[n] mu_stocks = theta_s[site] ./ omega_s[site];
  vector[n] mu_influx = theta_s[site];
  vector[n] mu_outflux = theta_s[site];
}
model {
  log(stocks) ~ normal(log(mu_stocks), sigma_stocks);
  log(influx) ~ normal(log(mu_influx), sigma_influx);
  log(outflux) ~ normal(log(mu_outflux), sigma_outflux);
  
  theta_s ~ normal(mu_theta, sigma_theta);
  omega_s ~ normal(mu_omega, sigma_omega);
  
  sigma_stocks ~ std_normal();
  sigma_influx ~ std_normal();
  sigma_outflux ~ std_normal();
  sigma_theta ~ std_normal();
  sigma_omega ~ std_normal();
}
