library(rethinking)
library(R.matlab)

fpath <- "Z:\\Jesus\\Jittering\\FULLYCurated\\09_190702_Jittering_3720_1520_1520 yes VPM- good\\KS2_newChanMap\\BI_tryout\\psth_control.mat"

dat <- R.matlab::readMat(fpath)

d <- list(
  #neuron = as.integer( dat$neuron ),
  nid = as.integer( 1:dat$Ncl ),
  
  idx1 = as.integer( ( ( 1:d$Ncl-1) * d$Nbins )+1 ),
  idx2 = as.integer( ( 1:d$Ncl )*d$Nbins ),
  #bin = as.integer( dat$bin ),
  counts = matrix( as.integer( dat$counts ) , ncol = ncol(dat$counts) ) ,
  Ncl = as.integer( dat$Ncl ),
  Ntrials = as.integer( dat$Ntrials ),
  Nbins = as.integer( dat$Nbins ),
  y = as.integer( dat$trig )
)

bayes_psth_psth <- "C:/Users/neuro/Documents/MATLAB/AuxiliaryFuncs/Stan/Bayes_PSTH.stan"

mPSTH.1 <- stan( bayes_psth_psth, cores=6, chains=4, data=d )

mPSTH.1 <- ulam(
  alist(
    y ~ binomial(1, p),
    logit(p) <- counts[ idx1[nid]:idx2[nid], 1:Ntrials] %*% alpha[,nid],
    
    transpars> matrix[40, nid]:alpha <-
      compose_noncentered( sigma_alpha, L_Rho_alpha, z_alpha ),
    matrix[nid, 40]:z_alpha ~ normal( 0 , 1 ),
    
    vector[40]:sigma_alpha ~ dexp(1),
    cholesky_factor_corr[40]:L_Rho_alpha ~ lkj_corr_cholesky( 2 ),
    
    gq> matrix[40,40]:Rho_alpha <<- Chol_to_Corr( L_Rho_alpha )
  ), data = d, chains = 4, cores = 4, log_lik = TRUE )

d$counts[1:d$Nbins,] %*%  ( 1:d$Ntrials )

m14.3 <- ulam(
  alist(
    L ~ binomial(1,p),
    logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id,tid],
    
    # adaptive priors - non-centered
    transpars> matrix[actor,4]:alpha <-
      compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
    transpars> matrix[block_id,4]:beta <-
      compose_noncentered( sigma_block , L_Rho_block , z_block ),
    matrix[4,actor]:z_actor ~ normal( 0 , 1 ),
    matrix[4,block_id]:z_block ~ normal( 0 , 1 ),
    
    # fixed priors
    g[tid] ~ normal(0,1),
    vector[4]:sigma_actor ~ dexp(1),
    cholesky_factor_corr[4]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
    vector[4]:sigma_block ~ dexp(1),
    cholesky_factor_corr[4]:L_Rho_block ~ lkj_corr_cholesky( 2 ),
    
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[4,4]:Rho_actor <<- Chol_to_Corr(L_Rho_actor),
    gq> matrix[4,4]:Rho_block <<- Chol_to_Corr(L_Rho_block)
  ) , data=dat , chains=4 , cores=4 , log_lik=TRUE )
