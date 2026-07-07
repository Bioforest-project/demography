functions {
  real stp_integral(real t0, real delta, real tau, real k) {
    real u0 = k * t0 / tau;
    real log_I = log(delta) + k + log(tau) - (k + 1) * log(k)
                 + lgamma(k + 1) + log(gamma_q(k + 1, u0));
    return exp(log_I);
  }
}
data {
  int<lower=0> n_equ; // number of equilibrium observations
  int<lower=0> n_log; // number of logged observations
  int<lower=0> s; //number of sites
  int<lower=0> p_log; //number of logged plots
  vector[n_equ] y_equ;
  vector[n_equ] in_equ;
  vector[n_equ] out_equ;
  vector[n_log] y_log;
  vector[n_log] in_log;
  vector[n_log] out_log;
  vector[n_log] time;  // recovery time
  vector[n_log] dt;
  vector[p_log] t0; // first year of post logging measurements
  vector[p_log] y0; // agb at first year of post logging measurements
  array[n_equ] int<lower=0, upper=s> site_equ;
  array[n_log] int<lower=0, upper=s> site_log;
  array[n_log] int<lower=0, upper=p_log> plot_log;
  array[p_log] int<lower=0, upper=s> site_plot;
  array[2] real lambda_bounds;
  array[2] real alpha_bounds;
  array[2] real delta_bounds;
  array[2] real tau_bounds;
  array[2] real k_bounds;
  array[2] real omega_bounds;
  array[2] real theta_bounds;
}
parameters {
  real<lower=lambda_bounds[1], upper=lambda_bounds[2]> mu_lambda; 
  real<lower=0> sigma_lambda;
  vector<lower=lambda_bounds[1], upper=lambda_bounds[2]>[s] lambda_s;
  real<lower=alpha_bounds[1], upper=alpha_bounds[2]> mu_alpha; 
  real<lower=0> sigma_alpha;
  vector<lower=alpha_bounds[1], upper=alpha_bounds[2]>[p_log] alpha_p;
  real<lower=delta_bounds[1], upper=delta_bounds[2]> mu_delta; 
  real<lower=0> sigma_delta;
  vector<lower=delta_bounds[1], upper=delta_bounds[2]>[p_log] delta_p;
  real<lower=tau_bounds[1], upper=tau_bounds[2]> mu_tau; 
  real<lower=0> sigma_tau;
  vector<lower=tau_bounds[1], upper=tau_bounds[2]>[s] tau_s;
  vector<lower=k_bounds[1], upper=k_bounds[2]>[s] k_s;
  real<lower=omega_bounds[1], upper=omega_bounds[2]> mu_omega; 
  real<lower=0> sigma_omega;
  vector<lower=omega_bounds[1], upper=omega_bounds[2]>[s] omega_s;
  real<lower=theta_bounds[1], upper=theta_bounds[2]> mu_theta; 
  real<lower=0> sigma_theta;
  vector<lower=theta_bounds[1], upper=theta_bounds[2]>[s] theta_s;
  real<lower=0> sigma_y;
  real<lower=0> sigma_in;
  real<lower=0> sigma_out;
  real<lower=0> sigma_diff;
}
transformed parameters {
  vector[n_log] mu_in_log = omega_s[site_log] .* y_log + 
    delta_p[plot_log] .* pow((time-dt./2) ./ tau_s[site_log], k_s[site_log]) .*
    exp(k_s[site_log] .* (1 - (time-dt./2) ./ tau_s[site_log]));
  vector[n_log] mu_out_log = omega_s[site_log] .* y_log + alpha_p[plot_log] .* 
    exp(-lambda_s[site_log].*(time-dt./2));
  vector[p_log] extra_prod;
  vector[p_log] diff_log;
  for (p in 1:p_log) {
    extra_prod[p] = stp_integral(t0[p], delta_p[p], tau_s[site_plot[p]], k_s[site_plot[p]]);
    diff_log[p] = extra_prod[p] -
      alpha_p[p] / lambda_s[site_plot[p]] * exp(-lambda_s[site_plot[p]] * t0[p]);
  }
}
model {
  log(in_log) ~ normal(log(mu_in_log), sigma_in);
  log(out_log) ~ normal(log(mu_out_log), sigma_out);
  diff_log ~ normal(theta_s[site_plot] - y0, sigma_diff);
  
  log(y_equ) ~ normal(log(theta_s[site_equ]), sigma_y);
  log(in_equ) ~ normal(log(theta_s[site_equ].*omega_s[site_equ]), sigma_in);
  log(out_equ) ~ normal(log(theta_s[site_equ].*omega_s[site_equ]), sigma_out);
  
  lambda_s ~ normal(mu_lambda, sigma_lambda) T[lambda_bounds[1], lambda_bounds[2]];
  alpha_p ~ normal(mu_alpha, sigma_alpha) T[alpha_bounds[1], alpha_bounds[2]];
  delta_p ~ normal(mu_delta, sigma_delta) T[delta_bounds[1], delta_bounds[2]];
  tau_s ~ normal(mu_tau, sigma_tau) T[tau_bounds[1], tau_bounds[2]];
  omega_s ~ normal(mu_omega, sigma_omega) T[omega_bounds[1], omega_bounds[2]];
  theta_s ~ normal(mu_theta, sigma_theta) T[theta_bounds[1], theta_bounds[2]];
  
  sigma_y ~ std_normal();
  sigma_in ~ std_normal();
  sigma_out ~ std_normal();
  sigma_diff ~ std_normal();
  sigma_lambda ~ std_normal();
  sigma_alpha ~ std_normal();
  sigma_delta ~ std_normal();
  sigma_tau ~ std_normal();
  sigma_omega ~ std_normal();
  sigma_theta ~ std_normal();
}
