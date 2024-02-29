//
// Ben Bales
//

data {
  int N; // Number of observations
  array[N] real y;
}
parameters {
  // Parameters of measurement model
  ordered[2] mu;
  array[2] real<lower=0> sigma;
  
  // Initial state
  simplex[2] rho;
  
  // Rows of the transition matrix
  simplex[2] t1;
  simplex[2] t2;
}
transformed parameters {
  matrix[2, 2] Gamma = rep_matrix(0, 2, 2);
  matrix[2, N] log_omega;
  
  // Build the transition matrix
  Gamma[1,  : ] = t1';
  Gamma[2,  : ] = t2';
  
  // Compute the log likelihoods in each possible state
  for (n in 1 : N) {
    // The observation model could change with n, or vary in a number of
    //  different ways (which is why log_omega is passed in as an argument)
    log_omega[1, n] = normal_lpdf(y[n] | mu[1], sigma[1] );
    log_omega[2, n] = normal_lpdf(y[n] | mu[2], sigma[2] );
  }
}
model {
  mu ~ normal(0, 2);
  sigma ~ exponential( 1 );
  
  rho ~ dirichlet( rep_vector(5, 2) );
  
  t1 ~ dirichlet( rep_vector(2, 2) );
  t2 ~ dirichlet( rep_vector(2, 2) );
  
  target += hmm_marginal(log_omega, Gamma, rho);
}
generated quantities {
  matrix[2, N] hidden_probs = hmm_hidden_state_prob(log_omega, Gamma, rho);
  array[N] int y_sim = hmm_latent_rng(log_omega, Gamma, rho);
}