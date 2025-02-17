functions {
  vector int_stocks(
    int n, 
    vector x, 
    vector tau, 
    vector lambda, 
    vector delta, 
    vector thetainf, 
    vector theta0, 
    vector omega, 
    vector y0
    ) {
      vector[n] k;
      vector[n] mu;
      k = y0 - thetainf ./ omega - ((theta0 - thetainf) .* lambda) ./ 
      (omega .* lambda - 1) - (16 * exp(2) .* delta .* thetainf) ./ 
      (tau .* tau .* (2 * omega - tau) .* (2 * omega - tau) .* 
      (2 * omega - tau));
      mu = thetainf ./ omega + k .* exp(-omega .* x) +
      ((theta0 - thetainf) .* lambda) ./ (omega .* lambda - 1) .* 
      exp(-lambda .* x) +  exp(2.0) .* exp(- x ./ tau) .* exp(- x ./ tau) .*
      ((2 * delta .* thetainf) ./ (tau .* tau .* (2 * omega - tau))) .*
      (x .* x - x * 4 ./ (2 * omega - tau) + 8 ./ ((2 * omega - tau) .*
      (2 * omega - tau)));
      return mu ;
    }
    
    vector productivity(
      vector x, 
      vector tau, 
      vector lambda, 
      vector delta, 
      vector thetainf, 
      vector theta0
      ) {
        return theta0 + (thetainf - theta0) .* (1 - exp(-lambda .* x)) + 
        thetainf .* delta .* ((x .* x) ./ (tau .* tau) .* exp(2 .* (1 - x ./ tau)));
      }
      
}
data {
  int<lower=0> n_rec; // obs reco
  int<lower=0> n_old; // obs old
  int<lower=0> n_pre; // obs prelog
  int<lower=0> n_site;
  int<lower=0> n_plot_rec;
  vector[n_rec] y_rec;
  vector[n_old] y_old;
  vector[n_pre] y_pre;
  vector[n_rec] in_rec;
  vector[n_old] in_old;
  vector[n_pre] in_pre;
  vector[n_rec] out_rec;
  vector[n_old] out_old;
  vector[n_pre] out_pre;
  vector[n_rec] time;
  array[n_rec] int<lower=0, upper=n_site> site_rec;
  array[n_old] int<lower=0, upper=n_site> site_old;
  array[n_pre] int<lower=0, upper=n_site> site_pre;
  array[n_rec] int<lower=0, upper=n_plot_rec> plot_rec;
  array[n_plot_rec] int<lower=0, upper=n_site> site_plot_rec;
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
  real<lower=0> sigma_y_old;
  real<lower=0> sigma_y_pre;
  real<lower=0> sigma_y_rec;
  real<lower=0> sigma_in_old;
  real<lower=0> sigma_in_pre;
  real<lower=0> sigma_in_rec;
  real<lower=0> sigma_out_old;
  real<lower=0> sigma_out_pre;
  real<lower=0> sigma_out_rec;
}
transformed parameters {
  vector[n_plot_rec] theta0_p = thetaInf_s[site_plot_rec] .* dist_p;
  
  vector[n_old] mu_y_old = thetaInf_s[site_old]./omega_s[site_old];
  vector[n_pre] mu_y_pre = thetaInf_s[site_pre]./omega_s[site_pre];
  
  vector[n_old] mu_in_old = thetaInf_s[site_old];
  vector[n_pre] mu_in_pre = thetaInf_s[site_pre];
  
  vector[n_old] mu_out_old = thetaInf_s[site_old];
  vector[n_pre] mu_out_pre = thetaInf_s[site_pre];
  
  vector[n_rec] mu_y_rec = int_stocks(n_rec, time, tau_s[site_rec], 
  lambda_p[plot_rec], delta_p[plot_rec], thetaInf_s[site_rec], 
  theta0_p[plot_rec], omega_s[site_rec], y0_p[plot_rec]);
  vector[n_rec] mu_in_rec = productivity(time, tau_s[site_rec], 
  lambda_p[plot_rec], delta_p[plot_rec], thetaInf_s[site_rec], 
  theta0_p[plot_rec]);
  vector[n_rec] mu_out_rec = omega_s[site_rec] .* int_stocks(n_rec, time, 
  tau_s[site_rec], lambda_p[plot_rec], delta_p[plot_rec], thetaInf_s[site_rec], 
  theta0_p[plot_rec], omega_s[site_rec], y0_p[plot_rec]);
}
model {
  log(y_old) ~ normal(log(mu_y_old), sigma_y_old);
  log(y_pre) ~ normal(log(mu_y_pre), sigma_y_pre);
  log(y_rec) ~ normal(log(mu_y_rec), sigma_y_rec);
  log(in_old) ~ normal(log(mu_in_old), sigma_in_old);
  log(in_pre) ~ normal(log(mu_in_pre), sigma_in_pre);
  log(in_rec) ~ normal(log(mu_in_rec), sigma_in_rec);
  log(out_old) ~ normal(log(mu_out_old), sigma_out_old);
  log(out_pre) ~ normal(log(mu_out_pre), sigma_out_pre);
  log(out_rec) ~ normal(log(mu_out_rec), sigma_out_rec);
  dist_p ~ cauchy(mu_dist, sigma_dist);
  thetaInf_s ~ cauchy(mu_thetaInf, sigma_thetaInf);
  lambda_p ~ cauchy(mu_lambda, sigma_lambda);
  delta_p ~ cauchy(mu_delta, sigma_delta);
  tau_s ~ cauchy(mu_tau, sigma_tau);
  omega_s ~ cauchy(mu_omega, sigma_omega);
  y0_p ~ cauchy(mu_y0, sigma_y0);
  sigma_y_old ~ std_normal();
  sigma_y_pre ~ std_normal();
  sigma_y_rec ~ std_normal();
  sigma_in_old ~ std_normal();
  sigma_in_pre ~ std_normal();
  sigma_in_rec ~ std_normal();
  sigma_out_old ~ std_normal();
  sigma_out_pre ~ std_normal();
  sigma_out_rec ~ std_normal();
  sigma_dist ~ std_normal();
  sigma_thetaInf ~ std_normal();
  sigma_lambda ~ std_normal();
  sigma_delta ~ std_normal();
  sigma_tau ~ std_normal();
  sigma_omega ~ std_normal();
  sigma_y0 ~ std_normal();
}
