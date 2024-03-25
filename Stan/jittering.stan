data{
    int N;
    int Nc;
    int Nb;
    int Nt;
    int Nn;
    array[N] int y;
    array[N] int trig;
    array[N] int bin_id;
    array[N] int neuron_id;
    array[N] int tid;
}
parameters{
     matrix[Nc,Nn] z_actor;
     matrix[Nc,Nb] z_block;
     matrix[Nc,Nt] z_eta;
     vector[Nc] g;
     vector<lower=0>[Nc] sigma_actor;
     cholesky_factor_corr[Nc] L_Rho_actor;
     vector<lower=0>[Nc] sigma_block;
     cholesky_factor_corr[Nc] L_Rho_block;
     vector<lower=0>[Nc] sigma_eta;
     cholesky_factor_corr[Nc] L_Rho_eta;
}
transformed parameters{
     matrix[Nn,Nc] alpha;
     matrix[Nb,Nc] beta;
     matrix[Nt,Nc] eta;
    eta = (diag_pre_multiply(sigma_eta, L_Rho_eta) * z_eta)';
    beta = (diag_pre_multiply(sigma_block, L_Rho_block) * z_block)';
    alpha = (diag_pre_multiply(sigma_actor, L_Rho_actor) * z_actor)';
}
model{
     vector[N] p;
     vector[N] lambda;
    L_Rho_eta ~ lkj_corr_cholesky( 2 );
    sigma_eta ~ exponential( 1 );
    L_Rho_block ~ lkj_corr_cholesky( 2 );
    sigma_block ~ exponential( 1 );
    L_Rho_actor ~ lkj_corr_cholesky( 2 );
    sigma_actor ~ exponential( 1 );
    g ~ normal( 0 , 1 );
    to_vector( z_eta ) ~ normal( 0 , 1 );
    to_vector( z_block ) ~ normal( 0 , 1 );
    to_vector( z_actor ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        lambda[i] = g[tid[i]] + alpha[neuron_id[i], tid[i]];
        lambda[i] = exp(lambda[i]);
    }
    for ( i in 1:N ) {
        p[i] = beta[bin_id[i], tid[i]] + eta[trig[i], tid[i]];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:N ) {
        if ( y[i]==0 )
            target += log_mix( p[i] , 0 , poisson_lpmf(0|lambda[i]) );
        if ( y[i] > 0 )
            target += log1m( p[i] ) + poisson_lpmf(y[i] | lambda[i] );
    }
}
generated quantities{
    vector[N] log_lik;
     vector[N] p;
     vector[N] lambda;
     matrix[Nc,Nc] Rho_actor;
     matrix[Nc,Nc] Rho_eta;
     matrix[Nc,Nc] Rho_block;
    Rho_block = multiply_lower_tri_self_transpose(L_Rho_block);
    Rho_eta = multiply_lower_tri_self_transpose(L_Rho_eta);
    Rho_actor = multiply_lower_tri_self_transpose(L_Rho_actor);
    for ( i in 1:N ) {
        lambda[i] = g[tid[i]] + alpha[neuron_id[i], tid[i]];
        lambda[i] = exp(lambda[i]);
    }
    for ( i in 1:N ) {
        p[i] = beta[bin_id[i], tid[i]] + eta[trig[i], tid[i]];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:N ) {
        if ( y[i]==0 )
            log_lik[i] = log_mix( p[i] , 0 , poisson_lpmf(0|lambda[i]) );
        if ( y[i] > 0 )
            log_lik[i] = log1m( p[i] ) + poisson_lpmf(y[i] | lambda[i] );
    }
}


