data {
  int<lower=1> n_equ;
  int<lower=1> n_log;
  int<lower=1> s;
  int<lower=1> j;
  int<lower=1> p_log;

  vector[n_equ] y_equ;
  vector[n_equ] in_equ;
  vector[n_equ] out_equ;

  vector[n_log] y_log;
  vector[n_log] in_log;
  vector[n_log] out_log;

  vector[n_log] time;
  vector[n_log] dt;

  array[n_equ] int<lower=1, upper=s> site_equ;
  array[n_log] int<lower=1, upper=s> site_log;
  array[n_log] int<lower=1, upper=p_log> plot_log;

  array[n_equ] int<lower=1, upper=j> size_equ;
  array[n_log] int<lower=1, upper=j> size_log;

  array[p_log] int<lower=1, upper=s> site_plot;

  array[2] real lambda_bounds;
  array[2] real alpha_bounds;
  array[2] real delta_bounds;
  array[2] real tau_bounds;
  array[2] real k_bounds;
  array[2] real omega_bounds;
  array[2] real theta_bounds;
}

parameters {
  // hierarchical means
  vector[j] mu_lambda;
  vector[j] mu_alpha;
  vector[j] mu_delta;
  vector[j] mu_tau;
  vector[j] mu_omega_p;
  vector[j-1] mu_omega_m;
  vector[j] mu_theta;

  // hierarchical SDs
  real<lower=0> sigma_lambda;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_delta;
  real<lower=0> sigma_tau;
  real<lower=0> sigma_omega;
  real<lower=0> sigma_theta;

  // varying effects
  matrix[s, j] lambda_s;
  matrix[p_log, j] alpha_p;
  matrix[p_log, j] delta_p;
  matrix[s, j] tau_s;
  matrix[s, j] omega_p_s;
  matrix[s, j-1] omega_m_s;
  matrix[s, j] theta_s;

  vector<lower=k_bounds[1], upper=k_bounds[2]>[j] k;

  // observation noise
  real<lower=0> sigma_y;
  real<lower=0> sigma_in;
  real<lower=0> sigma_out;
}

transformed parameters {
  vector[n_log] mu_in_log;
  vector[n_log] mu_out_log;
  vector[n_log] t_mid = time - dt * 0.5;
  vector[n_equ] theta_n; 
  vector[n_equ] omega_n; 
  
  matrix[s, j] omega_m_s_all;
  for (i in 1:(j-1)) {
    omega_m_s_all[,i] = omega_m_s[,i];
  }
  omega_m_s_all[,j] = sum()
  
  for (n in 1:n_equ) {
    theta_n[n] = theta_s[site_equ[n], size_equ[n]];
    omega_p_n[n] = omega_p_s[site_equ[n], size_equ[n]];
    omega_m_n[n] = omega_m_s_all[site_equ[n], size_equ[n]];
  }

  for (n in 1:n_log) {
    mu_in_log[n] =
      omega_s[site_log[n], size_log[n]] * y_log[n] +
      delta_p[plot_log[n], size_log[n]] *
      pow(t_mid[n] / tau_s[site_log[n], size_log[n]], k[size_log[n]]) *
      exp(k[size_log[n]] * (1 - t_mid[n] / tau_s[site_log[n], size_log[n]]));
    mu_out_log[n] =
      omega_s[site_log[n], size_log[n]] * y_log[n] +
      alpha_p[plot_log[n], size_log[n]] *
      exp(-lambda_s[site_log[n], size_log[n]] * t_mid[n]);
  }
}

model {
  // --- Likelihood (vectorized) ---
  log(in_log) ~ normal(log(mu_in_log), sigma_in);
  log(out_log) ~ normal(log(mu_out_log), sigma_out);

  log(y_equ) ~ normal(log(theta_n), sigma_y);
  log(in_equ) ~ normal(log(theta_n .* omega_p_n), sigma_in);
  log(out_equ) ~ normal(log(theta_n .* omega_m_n), sigma_out);

  // --- Hierarchical priors  ---
  for (jj in 1:j) {
    lambda_s[, jj] ~ normal(mu_lambda[jj], sigma_lambda);
    alpha_p[, jj]  ~ normal(mu_alpha[jj], sigma_alpha);
    delta_p[, jj]  ~ normal(mu_delta[jj], sigma_delta);
    tau_s[, jj]    ~ normal(mu_tau[jj], sigma_tau);
    omega_s[, jj]  ~ normal(mu_omega[jj], sigma_omega);
    theta_s[, jj]  ~ normal(mu_theta[jj], sigma_theta);
  }

  // --- Hyperpriors ---
  sigma_lambda ~ normal(0, 1);
  sigma_alpha  ~ normal(0, 1);
  sigma_delta  ~ normal(0, 1);
  sigma_tau    ~ normal(0, 1);
  sigma_omega  ~ normal(0, 1);
  sigma_theta  ~ normal(0, 1);
  sigma_y   ~ normal(0, 1);
  sigma_in  ~ normal(0, 1);
  sigma_out ~ normal(0, 1);
}
