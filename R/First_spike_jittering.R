library(rethinking)
library(R.matlab)

fpath <- "Z:\\Jesus\\Jittering\\FULLYCurated\\05_190701_Jesus_Emili_Jittering_3800_1600_1500 yes VPM good\\KS2_newChanMap\\FirstSpikes_perUnit_and_condition_filtered.mat"

dat <- R.matlab::readMat(fpath)

Z_fs <- standardize( dat$fs.puc[,3] )

d <- list(
  Fs = Z_fs,
  tid = as.integer( dat$fs.puc[,1] ),
  neuron = as.integer( dat$fs.puc[,2] )
  )

mFSJ.f <- ulam(
  alist(
    Fs ~ normal( mu, sigma ),
    # tid == 1 is the neuron with it's natural latency and trial-to-trial variance
    mu <- latency[neuron] + mod_lat[tid],
    log(sigma) <- jitt[neuron] + mod_jit[tid],
    
    #Adaptive, non-centered priors
    transpars> matrix[neuron,6]:latency <-
      compose_noncentered( sigma_latency , L_Rho_latency , z_latency ),
    matrix[6,neuron]:z_latency ~ normal( 0 , 1 ),
    
    transpars> matrix[neuron,6]:jitt <-
      compose_noncentered( sigma_jitt , L_Rho_jitt , z_jitt ),
    matrix[6,neuron]:z_jitt ~ normal( 0 , 1 ),
    
    #Fixed priors
    vector[6]:sigma_latency ~ exponential( 1 ),
    cholesky_factor_corr[6]:L_Rho_latency ~ lkj_corr_cholesky( 2 ),
    vector[6]:sigma_jitt ~ exponential( 2 ),
    cholesky_factor_corr[6]:L_Rho_jitt ~ lkj_corr_cholesky( 2 ),
    
    #Finally, compute ordinary correlation matrices from Cholesky factors
    gq> matrix[6,6]:Rho_latency <<- Chol_to_Corr( L_Rho_latency ),
    gq> matrix[6,6]:Rho_jitt <<- Chol_to_Corr( L_Rho_jitt )
  ), data = d , chains = 4, cores = 4)

mFSJ <- ulam(
  alist(
    Fs ~ normal( mu, sigma ),
    # tid == 1 is the neuron with it's natural latency and trial-to-trial variance
    mu <- latency[neuron,tid],
    log(sigma) <- jitt[neuron,tid],
    
    #Adaptive, non-centered priors
    transpars> matrix[neuron,6]:latency <-
      compose_noncentered( sigma_latency , L_Rho_latency , z_latency ),
    matrix[6,neuron]:z_latency ~ normal( 0 , 1 ),
    
    transpars> matrix[neuron,6]:jitt <-
      compose_noncentered( sigma_jitt , L_Rho_jitt , z_jitt ),
    matrix[6,neuron]:z_jitt ~ normal( 0 , 1 ),
    
    #Fixed priors
    vector[6]:sigma_latency ~ exponential( 1 ),
    cholesky_factor_corr[6]:L_Rho_latency ~ lkj_corr_cholesky( 2 ),
    vector[6]:sigma_jitt ~ exponential( 2 ),
    cholesky_factor_corr[6]:L_Rho_jitt ~ lkj_corr_cholesky( 2 ),
    
    #Finally, compute ordinary correlation matrices from Cholesky factors
    gq> matrix[6,6]:Rho_latency <<- Chol_to_Corr( L_Rho_latency ),
    gq> matrix[6,6]:Rho_jitt <<- Chol_to_Corr( L_Rho_jitt )
  ), data = d , chains = 4, cores = 4)

post <- link( mFSJ )
params <- extract.samples( mFSJ )

R.matlab::writeMat("Z:\\Jesus\\Jittering\\FULLYCurated\\05_190701_Jesus_Emili_Jittering_3800_1600_1500 yes VPM good\\KS2_newChanMap\\Jittering_BI.mat", post = post, params = params)