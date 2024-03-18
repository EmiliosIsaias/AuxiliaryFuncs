data{
     int Nbins;
     int Ntrials;
     int Ncl;
    array[Nbins] int y;
    array[Ncl] int nid;
    array[Ncl] int idx1;
    array[Ncl] int idx2;
     matrix[Nbins*Ncl,Ntrials] counts;
}
parameters{
     matrix[Ncl,Ntrials] z_alpha;
     vector<lower=0>[40] sigma_alpha;
     cholesky_factor_corr[40] L_Rho_alpha;
}
transformed parameters{
     matrix[Ntrials,Ncl] alpha;
    alpha = (diag_pre_multiply(sigma_alpha, L_Rho_alpha) * z_alpha)';
}
model{
     vector[Nbins, Ncl] p;
    L_Rho_alpha ~ lkj_corr_cholesky( 2 );
    sigma_alpha ~ exponential( 1 );
    to_vector( z_alpha ) ~ normal( 0 , 1 );
    for (i in 1:Ncl) {
        p[,i] = counts[idx1[nid[i]]:idx2[nid[i]], ] * alpha[, nid[i]] ;
        p[,i] = inv_logit(p[,i]);
    }
    y ~ binomial( 1 , p );
}
generated quantities{
    vector[Nbins] log_lik;
     vector[3546000] p;
     matrix[40,40] Rho_alpha;
    Rho_alpha = multiply_lower_tri_self_transpose(L_Rho_alpha);
    for ( i in 1:3546000 ) {
        p[i] = counts[idx1[nid[i]]:idx2[nid[i]], ] * alpha[, nid[i]];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:450 ) log_lik[i] = binomial_lpmf( y[i] | 1 , p[i] );
}


