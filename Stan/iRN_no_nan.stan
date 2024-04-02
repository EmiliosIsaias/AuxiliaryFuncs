data{
    int Ns;
    int Nc;
    int Nm;
    int N;
    array[N] y;
    array[N] int block_id;
    array[N] int actor_id;
    array[N] int tid;
}
parameters{
     matrix[Nc,Nm] z_actor;
     matrix[Nc,Ns] z_beta;
     vector[Nc] g;
     vector<lower=0>[Nc] sigma_actor;
     cholesky_factor_corr[Nc] L_Rho_actor;
     vector<lower=0>[Nc] sigma_beta;
     cholesky_factor_corr[Nc] L_Rho_beta;
     real<lower=0> sigma;
}
transformed parameters{
     matrix[Nm,Nc] alpha;
     matrix[Ns,Nc] beta;
    beta = (diag_pre_multiply(sigma_beta, L_Rho_beta) * z_beta)';
    alpha = (diag_pre_multiply(sigma_actor, L_Rho_actor) * z_actor)';
}
model{
     vector[N] mu;
    sigma ~ exponential( 1 );
    L_Rho_beta ~ lkj_corr_cholesky( 2 );
    sigma_beta ~ exponential( 1 );
    L_Rho_actor ~ lkj_corr_cholesky( 2 );
    sigma_actor ~ exponential( 1 );
    g ~ normal( 0 , 1 );
    to_vector( z_beta ) ~ normal( 0 , 1 );
    to_vector( z_actor ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        mu[i] = g[tid[i]] + alpha[actor_id[i], tid[i]] + beta[block_id[i], tid[i]];
    }
    y ~ normal( mu , sigma );
}
generated quantities{
    vector[N] log_lik;
     vector[N] mu;
     matrix[Nc,Nc] Rho_actor;
     matrix[Nc,Nc] Rho_beta;
    Rho_beta = multiply_lower_tri_self_transpose(L_Rho_beta);
    Rho_actor = multiply_lower_tri_self_transpose(L_Rho_actor);
    for ( i in 1:N ) {
        mu[i] = g[tid[i]] + alpha[actor_id[i], tid[i]] + beta[block_id[i], tid[i]];
    }
    for ( i in 1:N ) log_lik[i] = normal_lpdf( y[i] | mu[i] , sigma );
}