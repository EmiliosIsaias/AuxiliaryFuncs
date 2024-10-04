
data {
  int<lower=1> K; // K outcome dimensions
  int<lower=1> J; // J predictors
  int<lower=0> N; // N measurements in x
  array[N] vector[J] x; // NxJ matrix
  array[N] vector[K] y; // NxK matrix
}
parameters {
  matrix[K, J] beta;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] L_sigma;
  matrix[N, K] Z; // Standard normal variables
}

transformed parameters {
  array[N] vector[K] mu;
  matrix[K, K] L_Sigma;
  for (n in 1:N) {
    mu[n] = beta * x[n];
  }
  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
}

model {
  to_vector(beta) ~ normal(0, 2);
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ normal(0, 1);
  to_vector(Z) ~ normal(0, 1);
  for (n in 1:N) {
    y[n] ~ multi_normal_cholesky(mu[n], L_Sigma);
  }
}
