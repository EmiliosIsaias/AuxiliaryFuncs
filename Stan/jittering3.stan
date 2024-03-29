data{
    int N;
    int Nn;
    int Nb;
    int Nc;
    array[N] int trig;
    array[N] int psth;
    array[N] int bin_id;
    array[N] int neuron_id;
    array[N] int tid;
}
parameters{
     matrix[Nc,Nn] z_actor;
     matrix[Nc,Nb] z_block;
     vector[Nc] g;
     vector<lower=0>[Nc] sigma_actor;
     cholesky_factor_corr[Nc] L_Rho_actor;
     vector<lower=0>[Nc] sigma_block;
     cholesky_factor_corr[Nc] L_Rho_block;
}
transformed parameters{
     matrix[Nn,Nc] alpha;
     matrix[Nb,Nc] beta;
    beta = (diag_pre_multiply(sigma_block, L_Rho_block) * z_block)';
    alpha = (diag_pre_multiply(sigma_actor, L_Rho_actor) * z_actor)';
}
model{
     vector[N] p;
     vector[N] lambda;
    L_Rho_block ~ lkj_corr_cholesky( 2 );
    sigma_block ~ exponential( 1 );
    L_Rho_actor ~ lkj_corr_cholesky( 2 );
    sigma_actor ~ exponential( 1 );
    g ~ normal( 0 , 1 );
    to_vector( z_block ) ~ normal( 0 , 1 );
    to_vector( z_actor ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        lambda[i] = g[tid[i]] + alpha[neuron_id[i], tid[i]];
        lambda[i] = exp(lambda[i]);
    }
    for ( i in 1:N ) {
        p[i] = beta[bin_id[i], tid[i]];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:N ) {
        if ( psth[i]==0 )
            target += log_mix( p[i] , 0 , poisson_lpmf(0|lambda[i]) );
        if ( psth[i] > 0 )
            target += log1m( p[i] ) + poisson_lpmf(psth[i] | lambda[i] );
    }
}
generated quantities{
    vector[N] log_lik;
     vector[N] p;
     vector[N] lambda;
     matrix[Nc,Nc] Rho_actor;
     matrix[Nc,Nc] Rho_block;
    Rho_block = multiply_lower_tri_self_transpose(L_Rho_block);
    Rho_actor = multiply_lower_tri_self_transpose(L_Rho_actor);
    for ( i in 1:N ) {
        lambda[i] = g[tid[i]] + alpha[neuron_id[i], tid[i]];
        lambda[i] = exp(lambda[i]);
    }
    for ( i in 1:N ) {
        p[i] = beta[bin_id[i], tid[i]];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:N ) {
        if ( psth[i]==0 )
            log_lik[i] = log_mix( p[i] , 0 , poisson_lpmf(0|lambda[i]) );
        if ( psth[i] > 0 )
            log_lik[i] = log1m( p[i] ) + poisson_lpmf(psth[i] | lambda[i] );
    }
}


