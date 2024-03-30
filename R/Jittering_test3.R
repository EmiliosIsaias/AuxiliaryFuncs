library(rethinking)
library(R.matlab)

fpath <- "D:\\Linear_PSTH.mat"

dat <- R.matlab::readMat( fpath )

d <- list(
  bin_id = as.integer( dat$bin.id ),
  neuron_id = as.integer( dat$neuron.id ),
  tid = as.integer( dat$tid ),
  trig = as.integer( dat$y ),
  y = as.integer( dat$psth.linear )
)

mJitt.1 <- ulam(
  alist(
    y ~ dzipois(p, lambda),
    logit(p) <- g[tid] + alpha[neuron_id,tid] + beta[bin_id,tid] + eta[trig,tid],
    log(lambda) <- g[tid] + alpha[neuron_id,tid] + beta[bin_id,tid] + eta[trig,tid],
    
    # adaptive priors - non-centered
    transpars> matrix[neuron_id,6]:alpha <-
      compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
    transpars> matrix[bin_id,6]:beta <-
      compose_noncentered( sigma_block , L_Rho_block , z_block ),
    transpars> matrix[trig,6]:eta <-
      compose_noncentered( sigma_eta , L_Rho_eta , z_eta ),
    matrix[6,neuron_id]:z_actor ~ normal( 0 , 1 ),
    matrix[6,bin_id]:z_block ~ normal( 0 , 1 ),
    matrix[6,trig]:z_eta ~ normal( 0 , 1 ),
    
    # fixed priors
    g[tid] ~ normal( 0 , 1 ),
    vector[6]:sigma_actor ~ dexp(1),
    cholesky_factor_corr[6]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
    vector[6]:sigma_block ~ dexp(1),
    cholesky_factor_corr[6]:L_Rho_block ~ lkj_corr_cholesky( 2 ),
    vector[6]:sigma_eta ~ dexp(1),
    cholesky_factor_corr[6]:L_Rho_eta ~ lkj_corr_cholesky( 2 ),
    
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[6,6]:Rho_actor <<- Chol_to_Corr(L_Rho_actor),
    gq> matrix[6,6]:Rho_eta <<- Chol_to_Corr(L_Rho_eta),
    gq> matrix[6,6]:Rho_block <<- Chol_to_Corr(L_Rho_block)
  ), data = d, chains = 4, cores = 4, log_lik = TRUE )

mJitt.2 <- ulam(
  alist(
    y ~ dzipois(p, lambda),
    logit(p) <- beta[bin_id,tid] + eta[trig,tid],
    log(lambda) <- g[tid] + alpha[neuron_id,tid],
    
    # adaptive priors - non-centered
    transpars> matrix[neuron_id,6]:alpha <-
      compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
    transpars> matrix[bin_id,6]:beta <-
      compose_noncentered( sigma_block , L_Rho_block , z_block ),
    transpars> matrix[trig,6]:eta <-
      compose_noncentered( sigma_eta , L_Rho_eta , z_eta ),
    matrix[6,neuron_id]:z_actor ~ normal( 0 , 1 ),
    matrix[6,bin_id]:z_block ~ normal( 0 , 1 ),
    matrix[6,trig]:z_eta ~ normal( 0 , 1 ),
    
    # fixed priors
    g[tid] ~ normal( 0 , 1 ),
    vector[6]:sigma_actor ~ dexp(1),
    cholesky_factor_corr[6]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
    vector[6]:sigma_block ~ dexp(1),
    cholesky_factor_corr[6]:L_Rho_block ~ lkj_corr_cholesky( 2 ),
    vector[6]:sigma_eta ~ dexp(1),
    cholesky_factor_corr[6]:L_Rho_eta ~ lkj_corr_cholesky( 2 ),
    
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[6,6]:Rho_actor <<- Chol_to_Corr(L_Rho_actor),
    gq> matrix[6,6]:Rho_eta <<- Chol_to_Corr(L_Rho_eta),
    gq> matrix[6,6]:Rho_block <<- Chol_to_Corr(L_Rho_block)
  ), data = d, chains = 4, cores = 4, log_lik = TRUE )

mJitt.3 <- ulam(
  alist(
    
    y ~ dzipois(p, lambda),
    logit(p) <- beta[bin_id,tid],
    log(lambda) <- g[tid] + alpha[neuron_id,tid],
    
    # adaptive priors - non-centered
    transpars> matrix[neuron_id,6]:alpha <-
      compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
    transpars> matrix[bin_id,6]:beta <-
      compose_noncentered( sigma_block , L_Rho_block , z_block ),

    matrix[6,neuron_id]:z_actor ~ normal( 0 , 1 ),
    matrix[6,bin_id]:z_block ~ normal( 0 , 1 ),

    # fixed priors
    g[tid] ~ normal( 0 , 1 ),
    vector[6]:sigma_actor ~ dexp( 1 ),
    cholesky_factor_corr[6]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
    vector[6]:sigma_block ~ dexp( 1 ),
    cholesky_factor_corr[6]:L_Rho_block ~ lkj_corr_cholesky( 2 ),

    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[6,6]:Rho_actor <<- Chol_to_Corr(L_Rho_actor),
    gq> matrix[6,6]:Rho_block <<- Chol_to_Corr(L_Rho_block)
    
  ), data = d, chains = 2, cores = 4 , threads = 2)
