library(cmdstanr)
library(ggplot2)

N = 100 # 100 measurements
K = 3   # 3 states
states = rep(1, N)
states[1] = 1 # Start in state 1
for(n in 2:length(states)) {
  if(states[n - 1] == 1)
    states[n] = sample(c(1, 2), size = 1, prob = c(0.5, 0.5))
  else if(states[n - 1] == 2)
    states[n] = sample(c(1, 2, 3), size = 1, prob = c(0.25, 0.5, 0.25))
  else if(states[n - 1] == 3)
    states[n] = sample(c(2, 3), size = 1, prob = c(0.5, 0.5))
}

qplot(1:N, states)

mus = c(1.0, 5.0, 9.0)
sigma = 2.0
y = rnorm(N, mus[states], sd = sigma)

stan_model <- "data {
  int N; // Number of observations
  array[N] real y;
}
parameters {
  // Rows of the transition matrix
  simplex[2] t1;
  simplex[3] t2;
  simplex[2] t3;
  
  // Initial state
  simplex[3] rho;
  
  // Parameters of measurement model
  vector[3] mu;
  real<lower = 0.0> sigma;
}
transformed parameters {
  matrix[3, 3] gamma = rep_matrix(0, 3, 3);
  matrix[3, N] log_omega;
  
  // Build the transition matrix
  gamma[1, 1:2] = to_row_vector(t1);
  gamma[2, ] = to_row_vector(t2);
  gamma[3, 2:3] = to_row_vector(t3);
  
  // Compute the log likelihoods in each possible state
  for(n in 1:N) {
    // The observation model could change with n, or vary in a number of
    //  different ways (which is why log_omega is passed in as an argument)
    log_omega[1, n] = normal_lpdf(y[n] | mu[1], sigma);
    log_omega[2, n] = normal_lpdf(y[n] | mu[2], sigma);
    log_omega[3, n] = normal_lpdf(y[n] | mu[3], sigma);
  }
}
model {
  mu ~ normal(0, 1);
  sigma ~ normal(0, 1);

  target += hmm_marginal(log_omega, gamma, rho);
}"

model <- stan(model_code=stan_model, data=list(N=N, y=y), chains=3, cores=3)
