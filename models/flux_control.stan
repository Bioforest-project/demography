data {
  int<lower=0> n; // number of old-growth observations
  int<lower=0> s;
  vector[n] stocks;
  vector[n] influx;
  vector[n] outflux;
  array[n] int<lower=0, upper=s> site;   // site index
  array[2] real mu_thetaInf_bounds;
  array[2] real thetaInf_bounds;
  array[2] real omega_bounds;
}
parameters {
  real<lower=mu_thetaInf_bounds[1], upper=mu_thetaInf_bounds[2]> mu_thetaInf; 
  real<lower=0> sigma_thetaInf;
  vector<lower=thetaInf_bounds[1], upper=thetaInf_bounds[2]>[s] thetaInf_s;
  real<lower=omega_bounds[1], upper=omega_bounds[2]> mu_omega; // mortality rate
  real<lower=0> sigma_omega;
  vector<lower=omega_bounds[1], upper=omega_bounds[2]>[s] omega_s;
  real<lower=0> sigma_stocks;
  real<lower=0> sigma_influx;
  real<lower=0> sigma_outflux;
}
transformed parameters {
  vector[n] mu_stocks = thetaInf_s[site] ./ omega_s[site];
  vector[n] mu_influx = thetaInf_s[site];
  vector[n] mu_outflux = thetaInf_s[site];
}
model {
  log(stocks) ~ normal(log(mu_stocks), sigma_stocks);
  log(influx) ~ normal(log(mu_influx), sigma_influx);
  log(outflux) ~ normal(log(mu_outflux), sigma_outflux);
  
  thetaInf_s ~ normal(mu_thetaInf, sigma_thetaInf);
  omega_s ~ normal(mu_omega, sigma_omega);
  
  sigma_stocks ~ std_normal();
  sigma_influx ~ std_normal();
  sigma_outflux ~ std_normal();
  sigma_thetaInf ~ std_normal();
  sigma_omega ~ std_normal();
}
