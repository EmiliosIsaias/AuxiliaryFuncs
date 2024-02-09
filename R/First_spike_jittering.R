library(rethinking)
library(R.matlab)

fpath <- "Z:\\Jesus\\Jittering\\FULLYCurated\\05_190701_Jesus_Emili_Jittering_3800_1600_1500 yes VPM good\\KS2_newChanMap\\fs_puc_respFilt_tmFilt.mat"

dat <- R.matlab::readMat(fpath)

Z_fs <- standardize( dat$fs.puc[,3] )

d <- list(
  Fs = Z_fs,
  tid = as.integer( dat$fs.puc[,1] ),
  neuron = as.integer( dat$fs.puc[,2] )
  )

mFSJ.f1 <- ulam(
  alist(
    Fs ~ normal( mu, sigma ),
    
    mu <- g[tid] + alpha[neuron,tid],
    
    transpars> matrix[neuron,6]:alpha <- 
      compose_noncentered( sigma_alpha, L_Rho_alpha, z_alpha),
    matrix[6,neuron]:z_alpha ~ normal(0, 1),
    vector[6]:sigma_alpha ~ exponential( 1 ),
    cholesky_factor_corr[6]:L_Rho_alpha ~ lkj_corr_cholesky(3),
    
    g[tid] ~ normal(0, 1),
    sigma ~ exponential( 1 ),
    
    gq> matrix[6,6]:Rho_alpha <<- Chol_to_Corr(L_Rho_alpha)
  ), data = d , chains = 4, cores = 4)

mFSJ.f1_stan <- "data{
     vector[5912] Fs;
    array[5912] int neuron;
    array[5912] int tid;
}
parameters{
     matrix[6,164] z_alpha;
     vector<lower=0>[6] sigma_alpha;
     cholesky_factor_corr[6] L_Rho_alpha;
     vector[6] g;
     real<lower=0> sigma;
}
transformed parameters{
     matrix[164,6] alpha;
    alpha = (diag_pre_multiply(sigma_alpha, L_Rho_alpha) * z_alpha)';
}
model{
     vector[5912] mu;
    sigma ~ exponential( 1 );
    g ~ normal( 0 , 1 );
    L_Rho_alpha ~ lkj_corr_cholesky( 3 );
    sigma_alpha ~ exponential( 1 );
    to_vector( z_alpha ) ~ normal( 0 , 1 );
    for ( i in 1:5912 ) {
        mu[i] = g[tid[i]] + alpha[neuron[i], tid[i]];
    }
    Fs ~ normal( mu , sigma );
}
generated quantities{
     matrix[6,6] Rho_alpha;
    Rho_alpha = multiply_lower_tri_self_transpose(L_Rho_alpha);
}
"
mFSJ.f1 <- stan( model_code=mFSJ.f1_stan , data=d , chains=4 , cores=4 )

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

post <- link( mFSJ.f1 )
params <- extract.samples( mFSJ.f1 )

fs_att <- attributes(Z_fs)
fs_centre <- fs_att$`scaled:center`
fs_scale <- fs_att$`scaled:scale`

R.matlab::writeMat("C:\\Users\\jefe_\\mFSJ_f1_mu_respFilt_tmFilt.mat", post = post, params = params, fs_centre = fs_centre, fs_scale = fs_scale)
