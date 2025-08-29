data {
  int<lower=0> N;
  int<lower=0> K;  // number of 'humps' allowed
  int<lower=0> P;  
  vector[N] y;
  array[N] int<lower=0, upper=P> plot; 
  vector[N] time;
}
transformed data {
  vector[K] v_ones = rep_vector(1, K);
}
parameters {
  real<lower=0> sigma;
  vector<lower=log(5), upper=log(30)>[K] dlogtau;
  vector<lower=0,upper=5>[K] nu;
  matrix[P,K] delta;
  vector[K] deltaK;
  real<lower=0> sigma_delta;
}
transformed parameters {
  vector<lower=3,upper=50>[K] tau = cumulative_sum(exp(dlogtau));
  matrix[N,K] muk;
  for (k in 1:K) {
    for (n in 1:N) {
    muk[n,k] = delta[plot[n],k] * pow(time[n]/tau[k] * exp(1 - time[n]/tau[k]), exp(nu[k])) ;
    }
  }
  vector[N] mu = muk * v_ones;
}
model {
  y ~ normal(mu, sigma);
  
  //hyperdistribution of delta
  for (k in 1:K) {
    delta[,k] ~ normal(deltaK[k], sigma_delta);
  }
  
  dlogtau ~ normal(2, 2);
  nu ~ normal(2, 2);
  deltaK ~ normal(0, 5);
  sigma ~ normal(0, 10);
}

