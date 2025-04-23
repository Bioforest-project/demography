data {
  int<lower=0> n; // number of logged observations
  int<lower=0> s;
  int<lower=0> p;
  vector[n] stocks;
  vector[n] influx;
  vector[n] outflux;
  vector[n] time;  // recovery time
  array[n] int<lower=0, upper=s> site;   // site index
  array[n] int<lower=0, upper=p> plot;  // plot index 
  array[p] int<lower=0, upper=s> site_plot;
  vector[s] mu_thetaInf_s; // mean site asymptotic value in controls
  real<lower=0> sigma_thetaInf; // sd site asymptotic value in controls
  vector[s] mu_omega_s; // mean site asymptotic value in controls
  real<lower=0>  sigma_omega; // sd site asymptotic value in controls
  array[2] real thetaInf_bounds;
  array[2] real lambda_bounds;
  array[2] real dist_bounds;
  array[2] real delta_bounds;
  array[2] real tau_bounds;
  array[2] real omega_bounds;
  array[2] real y0_bounds;
}
parameters {
  real<lower=dist_bounds[1], upper=dist_bounds[2]> mu_dist; // starting point
  real<lower=0> sigma_dist;
  vector<lower=dist_bounds[1], upper=dist_bounds[2]>[p] dist_p;
  real<lower=lambda_bounds[1], upper=lambda_bounds[2]> mu_lambda; // rec rate
  real<lower=0> sigma_lambda;
  vector<lower=lambda_bounds[1], upper=lambda_bounds[2]>[p] lambda_p;
  vector<lower=thetaInf_bounds[1], upper=thetaInf_bounds[2]>[s] thetaInf_s;
  real<lower=delta_bounds[1], upper=delta_bounds[2]> mu_delta; // str var
  real<lower=0> sigma_delta;
  vector<lower=delta_bounds[1], upper=delta_bounds[2]>[p] delta_p;
  real<lower=tau_bounds[1], upper=tau_bounds[2]> mu_tau; //str time
  real<lower=0> sigma_tau;
  vector<lower=tau_bounds[1], upper=tau_bounds[2]>[s] tau_s;
  vector<lower=omega_bounds[1], upper=omega_bounds[2]>[s] omega_s;
  real<lower=y0_bounds[1], upper=y0_bounds[2]> mu_y0; // post-logging stocks
  real<lower=0> sigma_y0;
  vector<lower=y0_bounds[1], upper=y0_bounds[2]>[p] y0_p;
  real<lower=0> sigma_stocks;
  real<lower=0> sigma_influx;
  real<lower=0> sigma_outflux;
}
transformed parameters {
  vector[p] theta0_p = thetaInf_s[site_plot] .* dist_p;
  
  vector[n] mu_stocks = thetaInf_s[site] ./ omega_s[site] + 
  (y0_p[plot] - thetaInf_s[site] ./ omega_s[site] - 
  ((theta0_p[plot] - thetaInf_s[site]) .* lambda_p[plot]) ./ 
  (omega_s[site] .* lambda_p[plot] - 1) - 
  (16 * exp(2) .* delta_p[plot] .* thetaInf_s[site]) ./ 
  (tau_s[site] .* tau_s[site] .* 
  (2 * omega_s[site] - tau_s[site]) .* 
  (2 * omega_s[site] - tau_s[site]) .* 
  (2 * omega_s[site] - tau_s[site]))) .* 
  exp(-omega_s[site] .* time) +
  ((theta0_p[plot] - thetaInf_s[site]) .* lambda_p[plot]) ./ 
  (omega_s[site] .* lambda_p[plot] - 1) .* 
  exp(-lambda_p[plot] .* time) +  exp(2) .* 
  exp(- time ./ tau_s[site]) .* 
  exp(- time ./ tau_s[site]) .*
  ((2 * delta_p[plot] .* thetaInf_s[site]) ./
  (tau_s[site] .* tau_s[site] .* 
  (2 * omega_s[site] - tau_s[site]))) .*
  (time .* time - time * 4 ./ (2 * omega_s[site] - tau_s[site]) + 8 ./ 
  ((2 * omega_s[site] - tau_s[site]) .*
  (2 * omega_s[site] - tau_s[site])));
    
  vector[n] mu_influx = thetaInf_s[site] + 
    (theta0_p[plot] - thetaInf_s[site]) .* 
    exp(-lambda_p[plot] .* time) + 
    thetaInf_s[site] .* delta_p[plot] .* ((time .* time) ./ 
    (tau_s[site] .* tau_s[site]) .*
    exp(2 .* (1 - time ./ tau_s[site])));
  
  vector[n] mu_outflux = omega_s[site] .* (thetaInf_s[site] ./ 
  omega_s[site] + 
  (y0_p[plot] - thetaInf_s[site] ./ omega_s[site] - 
  ((theta0_p[plot] - thetaInf_s[site]) .* lambda_p[plot]) ./ 
  (omega_s[site] .* lambda_p[plot] - 1) - 
  (16 * exp(2) .* delta_p[plot] .* thetaInf_s[site]) ./ 
  (tau_s[site] .* tau_s[site] .* 
  (2 * omega_s[site] - tau_s[site]) .* 
  (2 * omega_s[site] - tau_s[site]) .* 
  (2 * omega_s[site] - tau_s[site]))) .* 
  exp(-omega_s[site] .* time) +
  ((theta0_p[plot] - thetaInf_s[site]) .* lambda_p[plot]) ./ 
  (omega_s[site] .* lambda_p[plot] - 1) .* 
  exp(-lambda_p[plot] .* time) +  exp(2) .* 
  exp(- time ./ tau_s[site]) .* 
  exp(- time ./ tau_s[site]) .* 
  ((2 * delta_p[plot] .* thetaInf_s[site]) ./ 
  (tau_s[site] .* tau_s[site] .* 
  (2 * omega_s[site] - tau_s[site]))) .*
  (time .* time - time * 4 ./ (2 * omega_s[site] - tau_s[site]) + 8 ./ 
  ((2 * omega_s[site] - tau_s[site]) .*
  (2 * omega_s[site] - tau_s[site])))
  );
}
model {
  log(stocks) ~ normal(log(mu_stocks), sigma_stocks);
  log(influx) ~ normal(log(mu_influx), sigma_influx);
  log(outflux) ~ normal(log(mu_outflux), sigma_outflux);
  
  // prior informed by control plots
  for(i in 1:s){
    thetaInf_s[i] ~ normal(mu_thetaInf_s[i], sigma_thetaInf);
    omega_s[i] ~ normal(mu_omega_s[i], sigma_omega);
  }
  
  dist_p ~ normal(mu_dist, sigma_dist);
  lambda_p ~ normal(mu_lambda, sigma_lambda);
  delta_p ~ normal(mu_delta, sigma_delta);
  tau_s ~ normal(mu_tau, sigma_tau);
  y0_p ~ normal(mu_y0, sigma_y0);
  
  sigma_stocks ~ std_normal();
  sigma_influx ~ std_normal();
  sigma_outflux ~ std_normal();
  sigma_dist ~ std_normal();
  sigma_lambda ~ std_normal();
  sigma_delta ~ std_normal();
  sigma_tau ~ std_normal();
  sigma_omega ~ std_normal();
  sigma_y0 ~ std_normal();
}
