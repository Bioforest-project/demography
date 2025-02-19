data {
  int<lower=0> n_rec; // number of logged observations
  int<lower=0> n_old; // number of old-growth observations
  int<lower=0> n_site;
  int<lower=0> n_plot_rec;
  vector[n_rec] stocks_rec;
  vector[n_rec] influx_rec;
  vector[n_rec] outflux_rec;
  vector[n_old] stocks_old;
  vector[n_old] influx_old;
  vector[n_old] outflux_old;
  vector[n_rec] time;  // recovery time
  array[n_old] int<lower=0, upper=n_site> site_old;   // site index
  array[n_rec] int<lower=0, upper=n_site> site_rec;   // site index
  array[n_rec] int<lower=0, upper=n_plot_rec> plot_rec;  // plot index 
  array[n_plot_rec] int<lower=0, upper=n_site> site_plot;
  array[2] real mu_thetaInf_bounds;
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
  vector<lower=dist_bounds[1], upper=dist_bounds[2]>[n_plot_rec] dist_p;
  real<lower=lambda_bounds[1], upper=lambda_bounds[2]> mu_lambda; // rec rate
  real<lower=0> sigma_lambda;
  vector<lower=lambda_bounds[1], upper=lambda_bounds[2]>[n_plot_rec] lambda_p;
  real<lower=mu_thetaInf_bounds[1], upper=mu_thetaInf_bounds[2]> mu_thetaInf; 
  real<lower=0> sigma_thetaInf;
  vector<lower=thetaInf_bounds[1], upper=thetaInf_bounds[2]>[n_site] thetaInf_s;
  real<lower=delta_bounds[1], upper=delta_bounds[2]> mu_delta; // str var
  real<lower=0> sigma_delta;
  vector<lower=delta_bounds[1], upper=delta_bounds[2]>[n_plot_rec] delta_p;
  real<lower=tau_bounds[1], upper=tau_bounds[2]> mu_tau; //str time
  real<lower=0> sigma_tau;
  vector<lower=tau_bounds[1], upper=tau_bounds[2]>[n_site] tau_s;
  real<lower=omega_bounds[1], upper=omega_bounds[2]> mu_omega; // mortality rate
  real<lower=0> sigma_omega;
  vector<lower=omega_bounds[1], upper=omega_bounds[2]>[n_site] omega_s;
  real<lower=y0_bounds[1], upper=y0_bounds[2]> mu_y0; // post-logging stocks
  real<lower=0> sigma_y0;
  vector<lower=y0_bounds[1], upper=y0_bounds[2]>[n_plot_rec] y0_p;
  real<lower=0> sigma_stocks;
  real<lower=0> sigma_influx;
  real<lower=0> sigma_outflux;
}
transformed parameters {
  vector[n_plot_rec] theta0_p = thetaInf_s[site_plot] .* dist_p;
  
  vector[n_old] mu_stocks_old = thetaInf_s[site_old] ./ omega_s[site_old];
  vector[n_old] mu_influx_old = thetaInf_s[site_old];
  vector[n_old] mu_outflux_old = thetaInf_s[site_old];
  
  vector[n_rec] mu_stocks_rec = thetaInf_s[site_rec] ./ omega_s[site_rec] + 
  (y0_p[plot_rec] - thetaInf_s[site_rec] ./ omega_s[site_rec] - 
  ((theta0_p[plot_rec] - thetaInf_s[site_rec]) .* lambda_p[plot_rec]) ./ 
  (omega_s[site_rec] .* lambda_p[plot_rec] - 1) - 
  (16 * exp(2) .* delta_p[plot_rec] .* thetaInf_s[site_rec]) ./ 
  (tau_s[site_rec] .* tau_s[site_rec] .* 
  (2 * omega_s[site_rec] - tau_s[site_rec]) .* 
  (2 * omega_s[site_rec] - tau_s[site_rec]) .* 
  (2 * omega_s[site_rec] - tau_s[site_rec]))) .* 
  exp(-omega_s[site_rec] .* time) +
  ((theta0_p[plot_rec] - thetaInf_s[site_rec]) .* lambda_p[plot_rec]) ./ 
  (omega_s[site_rec] .* lambda_p[plot_rec] - 1) .* 
  exp(-lambda_p[plot_rec] .* time) +  exp(2) .* 
  exp(- time ./ tau_s[site_rec]) .* 
  exp(- time ./ tau_s[site_rec]) .*
  ((2 * delta_p[plot_rec] .* thetaInf_s[site_rec]) ./
  (tau_s[site_rec] .* tau_s[site_rec] .* 
  (2 * omega_s[site_rec] - tau_s[site_rec]))) .*
  (time .* time - time * 4 ./ (2 * omega_s[site_rec] - tau_s[site_rec]) + 8 ./ 
  ((2 * omega_s[site_rec] - tau_s[site_rec]) .*
  (2 * omega_s[site_rec] - tau_s[site_rec])));
    
    vector[n_rec] mu_influx_rec = thetaInf_s[site_rec] + 
    (theta0_p[plot_rec] - thetaInf_s[site_rec]) .* 
    exp(-lambda_p[plot_rec] .* time) + 
    thetaInf_s[site_rec] .* delta_p[plot_rec] .* ((time .* time) ./ 
    (tau_s[site_rec] .* tau_s[site_rec]) .*
    exp(2 .* (1 - time ./ tau_s[site_rec])));
  
  vector[n_rec] mu_outflux_rec = omega_s[site_rec] .* (thetaInf_s[site_rec] ./ 
  omega_s[site_rec] + 
  (y0_p[plot_rec] - thetaInf_s[site_rec] ./ omega_s[site_rec] - 
  ((theta0_p[plot_rec] - thetaInf_s[site_rec]) .* lambda_p[plot_rec]) ./ 
  (omega_s[site_rec] .* lambda_p[plot_rec] - 1) - 
  (16 * exp(2) .* delta_p[plot_rec] .* thetaInf_s[site_rec]) ./ 
  (tau_s[site_rec] .* tau_s[site_rec] .* 
  (2 * omega_s[site_rec] - tau_s[site_rec]) .* 
  (2 * omega_s[site_rec] - tau_s[site_rec]) .* 
  (2 * omega_s[site_rec] - tau_s[site_rec]))) .* 
  exp(-omega_s[site_rec] .* time) +
  ((theta0_p[plot_rec] - thetaInf_s[site_rec]) .* lambda_p[plot_rec]) ./ 
  (omega_s[site_rec] .* lambda_p[plot_rec] - 1) .* 
  exp(-lambda_p[plot_rec] .* time) +  exp(2) .* 
  exp(- time ./ tau_s[site_rec]) .* 
  exp(- time ./ tau_s[site_rec]) .* 
  ((2 * delta_p[plot_rec] .* thetaInf_s[site_rec]) ./ 
  (tau_s[site_rec] .* tau_s[site_rec] .* 
  (2 * omega_s[site_rec] - tau_s[site_rec]))) .*
  (time .* time - time * 4 ./ (2 * omega_s[site_rec] - tau_s[site_rec]) + 8 ./ 
  ((2 * omega_s[site_rec] - tau_s[site_rec]) .*
  (2 * omega_s[site_rec] - tau_s[site_rec])))
  );
}
model {
  log(stocks_old) ~ normal(log(mu_stocks_old), sigma_stocks);
  log(influx_old) ~ normal(log(mu_influx_old), sigma_influx);
  log(outflux_old) ~ normal(log(mu_outflux_old), sigma_outflux);
  log(stocks_rec) ~ normal(log(mu_stocks_rec), sigma_stocks);
  log(influx_rec) ~ normal(log(mu_influx_rec), sigma_influx);
  log(outflux_rec) ~ normal(log(mu_outflux_rec), sigma_outflux);
  
  dist_p ~ cauchy(mu_dist, sigma_dist);
  thetaInf_s ~ cauchy(mu_thetaInf, sigma_thetaInf);
  lambda_p ~ cauchy(mu_lambda, sigma_lambda);
  delta_p ~ cauchy(mu_delta, sigma_delta);
  tau_s ~ cauchy(mu_tau, sigma_tau);
  omega_s ~ cauchy(mu_omega, sigma_omega);
  y0_p ~ cauchy(mu_y0, sigma_y0);
  
  sigma_stocks ~ std_normal();
  sigma_influx ~ std_normal();
  sigma_outflux ~ std_normal();
  sigma_dist ~ std_normal();
  sigma_thetaInf ~ std_normal();
  sigma_lambda ~ std_normal();
  sigma_delta ~ std_normal();
  sigma_tau ~ std_normal();
  sigma_omega ~ std_normal();
  sigma_y0 ~ std_normal();
}
