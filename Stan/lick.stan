data{
     int Nts;
     int Nm;
     int N;
    array[N] int y;
    array[N] int progress;
    array[N] int block_id;
    array[N] int actor;
    array[N] int tid;
}
parameters{
     matrix[Nts,Nm] z_actor;
     matrix[Nts,Nts] z_beta;
     matrix[Nts,1000] z_eta;
     vector[Nts] g;
     vector<lower=0>[Nts] sigma_actor;
     cholesky_factor_corr[Nts] L_Rho_actor;
     vector<lower=0>[Nts] sigma_beta;
     cholesky_factor_corr[Nts] L_Rho_beta;
     vector<lower=0>[Nts] sigma_eta;
     cholesky_factor_corr[Nts] L_Rho_eta;
     real<lower=0> sigma;
}
transformed parameters{
     matrix[Nm,Nts] alpha;
     matrix[Nts,Nts] beta;
     matrix[1000,Nts] eta;
    eta = (diag_pre_multiply(sigma_eta, L_Rho_eta) * z_eta)';
    beta = (diag_pre_multiply(sigma_beta, L_Rho_beta) * z_beta)';
    alpha = (diag_pre_multiply(sigma_actor, L_Rho_actor) * z_actor)';
}
model{
     vector[N] p;
    sigma ~ exponential( 1 );
    L_Rho_eta ~ lkj_corr_cholesky( 2 );
    sigma_eta ~ exponential( 1 );
    L_Rho_beta ~ lkj_corr_cholesky( 2 );
    sigma_beta ~ exponential( 1 );
    L_Rho_actor ~ lkj_corr_cholesky( 2 );
    sigma_actor ~ exponential( 1 );
    g ~ normal( 0 , 1 );
    to_vector( z_eta ) ~ normal( 0 , 1 );
    to_vector( z_beta ) ~ normal( 0 , 1 );
    to_vector( z_actor ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        p[i] = g[tid[i]] + alpha[actor[i], tid[i]] + beta[block_id[i], tid[i]] + eta[progress[i], tid[i]];
        p[i] = inv_logit(p[i]);
    }
    y ~ binomial( 1 , p );
}
generated quantities{
    vector[N] log_lik;
     vector[N] p;
     matrix[Nts,Nts] Rho_actor;
     matrix[Nts,Nts] Rho_beta;
     matrix[Nts,Nts] Rho_eta;
    Rho_eta = multiply_lower_tri_self_transpose(L_Rho_eta);
    Rho_beta = multiply_lower_tri_self_transpose(L_Rho_beta);
    Rho_actor = multiply_lower_tri_self_transpose(L_Rho_actor);
    for ( i in 1:N ) {
        p[i] = g[tid[i]] + alpha[actor[i], tid[i]] + beta[block_id[i], tid[i]] + eta[progress[i], tid[i]];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:N ) log_lik[i] = binomial_lpmf( y[i] | 1 , p[i] );
}


