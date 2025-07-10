data {
  int<lower=0> n_equ; // number of equilibrium observations
  int<lower=0> n_log; // number of logged observations
  int<lower=0> s; //number of sites
  // int<lower=0> p; //number of equilibrium plots
  int<lower=0> p_log; //number of logged plots
  vector[n_equ] y_equ;
  vector[n_equ] in_equ;
  vector[n_equ] out_equ;
  vector[n_log] y_log;
  vector[n_log] in_log;
  vector[n_log] out_log;
  vector[n_log] time;  // recovery time
  array[n_equ] int<lower=0, upper=s> site_equ;   // equilibrium site index
  array[n_log] int<lower=0, upper=s> site_log;   // logged site index
  array[n_log] int<lower=0, upper=p_log> plot_log;  // logged plot index 
  array[p_log] int<lower=0, upper=s> site_plot;  // logged plot to site index 
  array[2] real y0_bounds;
  array[2] real lambda_bounds;
  array[2] real alpha_bounds;
  array[2] real delta_bounds;
  array[2] real tau_bounds;
  array[2] real omega_bounds;
  array[2] real theta_bounds;
}
parameters {
  // post-logging stocks (logged plot)
  real<lower=y0_bounds[1], upper=y0_bounds[2]> mu_y0; 
  real<lower=0> sigma_y0;
  vector<lower=y0_bounds[1], upper=y0_bounds[2]>[p_log] y0_p;
  // mortality decrease rate (site)
  real<lower=lambda_bounds[1], upper=lambda_bounds[2]> mu_lambda; 
  real<lower=0> sigma_lambda;
  vector<lower=lambda_bounds[1], upper=lambda_bounds[2]>[s] lambda_s;
  // initial overmortality (logged plot)
  real<lower=alpha_bounds[1], upper=alpha_bounds[2]> mu_alpha; 
  real<lower=0> sigma_alpha;
  vector<lower=alpha_bounds[1], upper=alpha_bounds[2]>[p_log] alpha_p;
  // productivity stp height (logged plot)
  real<lower=delta_bounds[1], upper=delta_bounds[2]> mu_delta; 
  real<lower=0> sigma_delta;
  vector<lower=delta_bounds[1], upper=delta_bounds[2]>[p_log] delta_p;
  // productivity stp timing (site)
  real<lower=tau_bounds[1], upper=tau_bounds[2]> mu_tau; 
  real<lower=0> sigma_tau;
  vector<lower=tau_bounds[1], upper=tau_bounds[2]>[s] tau_s;
  // flux rate (site)
  real<lower=omega_bounds[1], upper=omega_bounds[2]> mu_omega; 
  real<lower=0> sigma_omega;
  vector<lower=omega_bounds[1], upper=omega_bounds[2]>[s] omega_s;
  // equilibrium stocks (site)
  real<lower=theta_bounds[1], upper=theta_bounds[2]> mu_theta; 
  real<lower=0> sigma_theta;
  vector<lower=theta_bounds[1], upper=theta_bounds[2]>[s] theta_s;
  // standard deviations
  real<lower=0> sigma_y;
  real<lower=0> sigma_in;
  real<lower=0> sigma_out;
}
transformed parameters {
  vector[p_log] theta_log = y0_p + exp(2) ./ 4 .* delta_p .* tau_s[site_plot] - 
  alpha_p ./ lambda_s[site_plot];
  vector[n_log] mu_y_log = y0_p[plot_log]+delta_p[plot_log] .* tau_s[site_log].* 
  exp(2) ./ 4 - delta_p[plot_log] .* exp((1 - time ./ tau_s[site_log]) .* 2) .*
  (time.*time ./ (2 .* tau_s[site_log]) + time ./ 2 + tau_s[site_log] ./ 4) -
  alpha_p[plot_log]./lambda_s[site_log].*(1-exp(-lambda_s[site_log].*time));
  vector[n_log] mu_in_log = omega_s[site_log] .* mu_y_log + 
  delta_p[plot_log] .*((time .* time) ./ (tau_s[site_log].*tau_s[site_log]) .*
  exp(2 .* (1 - time ./ tau_s[site_log])));
  vector[n_log] mu_out_log = omega_s[site_log] .* mu_y_log + alpha_p[plot_log] .* 
  exp(1-lambda_s[site_log].*time);
}
model {
  log(y_log) ~ normal(log(mu_y_log), sigma_y);
  log(in_log) ~ normal(log(mu_in_log), sigma_in);
  log(out_log) ~ normal(log(mu_out_log), sigma_out);
  
  log(theta_log) ~ normal(log(theta_s[site_plot]), sigma_theta);
  
  log(y_equ) ~ normal(log(theta_s[site_equ]), sigma_y);
  log(in_equ) ~ normal(log(theta_s[site_equ].*omega_s[site_log]), sigma_in);
  log(out_equ) ~ normal(log(theta_s[site_equ].*omega_s[site_log]), sigma_out);
  
  lambda_s ~ normal(mu_lambda, sigma_lambda);
  alpha_p ~ normal(mu_alpha, sigma_alpha);
  delta_p ~ normal(mu_delta, sigma_delta);
  tau_s ~ normal(mu_tau, sigma_tau);
  omega_s ~ normal(mu_omega, sigma_omega);
  y0_p ~ normal(mu_y0, sigma_y0);
  
  sigma_y ~ std_normal();
  sigma_in ~ std_normal();
  sigma_out ~ std_normal();
  sigma_lambda ~ std_normal();
  sigma_alpha ~ std_normal();
  sigma_delta ~ std_normal();
  sigma_tau ~ std_normal();
  sigma_omega ~ std_normal();
  sigma_y0 ~ std_normal();
}
