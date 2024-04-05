library(rethinking)
library(R.matlab)

fpath <- file.path("C:", "Users", "jefe_", "seadrive_root", "Emilio U", 
                   "Für meine Gruppen", "GDrive GrohLab", "Projects", 
                   "00 Salience", "Bayes Model Data", "Lick_data_4_R.mat", 
                   fsep = .Platform$file.sep )
dat <- R.matlab::readMat(fpath)

lick_mu <- mean( dat$lick.lat, na.rm = TRUE)
lick_sig <- sd( dat$lick.lat, na.rm = TRUE)

z_lick <- ( dat$lick.lat - lick_mu ) / lick_sig

d <- list(
  y = as.integer( dat$lick.flag ),
  #l_lat = z_lick,
  actor = as.integer( dat$mouse.id ),
  tid = as.integer( dat$tid ),
  block_id = as.integer( dat$block.id ),
  progress = as.integer( dat$progress.id ),
  N = as.integer( dat$N ),
  Nm = as.integer( dat$Nm ),
  Nts = as.integer( dat$Nts )
)

# d <- dat$mice.habituation[complete.cases(dat$mice.habituation),]
# B_z <- standardize(d[,4])
# d <- list(
#   actor = as.integer(d[,1]),
#   block_id = as.integer(d[,2]),
#   tid = as.integer(2L + d[,3]*2L),
#   B = B_z
# )
# d$tid[d$tid == 0L] = 1L

mLick.2 <- ulam( 
  alist(
    y ~ binomial( 1 , p ),
    logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id,tid],
    
    #Adaptive priors
    transpars> matrix[actor,4]:alpha <-
      compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
    transpars> matrix[block_id,4]:beta <-
      compose_noncentered( sigma_beta , L_Rho_beta , z_beta ),
    matrix[4,actor]:z_actor ~ normal( 0 , 1 ),
    matrix[4,block_id]:z_beta ~ normal( 0 , 1 ),

    #Fixed priors
    g[tid] ~ normal( 0 , 1 ),
    vector[4]:sigma_actor ~ exponential( 1 ),
    cholesky_factor_corr[4]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
    vector[4]:sigma_beta ~ exponential( 1 ),
    cholesky_factor_corr[4]:L_Rho_beta ~ lkj_corr_cholesky( 2 ),
    sigma ~ exponential( 1 ),

    #Finally, compute ordinary correlation matrices from Cholesky factors
    gq> matrix[4,4]:Rho_actor <<- Chol_to_Corr( L_Rho_actor ),
    gq> matrix[4,4]:Rho_beta <<- Chol_to_Corr( L_Rho_beta )
  ), data = d, chains = 4, cores = 4, log_lik = TRUE )

post <- extract.samples( mLick.2 )
pred <- link( mLick.2 )

fpath_out <- "D:\\lick_post_mLick2.mat"
R.matlab::writeMat(fpath_out, post = post)

fpath_out <- "D:\\lick_pred_mLick2.mat"
R.matlab::writeMat(fpath_out, pred = pred)


mLick.3 <- ulam( 
  alist(
    y ~ binomial( 1 , p ),
    logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id,tid] + 
      eta[progress, tid],
    
    #Adaptive priors
    transpars> matrix[actor,Nts]:alpha <-
      compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
    transpars> matrix[block_id,Nts]:beta <-
      compose_noncentered( sigma_beta , L_Rho_beta , z_beta ),
    transpars> matrix[progress,Nts]:eta <-
      compose_noncentered( sigma_eta , L_Rho_eta , z_eta ),
    matrix[Nts,actor]:z_actor ~ normal( 0 , 1 ),
    matrix[Nts,block_id]:z_beta ~ normal( 0 , 1 ),
    matrix[Nts,progress]:z_eta ~ normal( 0 , 1 ),
    
    #Fixed priors
    g[tid] ~ normal( 0 , 1 ),
    vector[Nts]:sigma_actor ~ exponential( 1 ),
    cholesky_factor_corr[Nts]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
    vector[Nts]:sigma_beta ~ exponential( 1 ),
    cholesky_factor_corr[Nts]:L_Rho_beta ~ lkj_corr_cholesky( 2 ),
    vector[Nts]:sigma_eta ~ exponential( 1 ),
    cholesky_factor_corr[Nts]:L_Rho_eta ~ lkj_corr_cholesky( 2 ),
    
    sigma ~ exponential( 1 ),
    
    #Finally, compute ordinary correlation matrices from Cholesky factors
    gq> matrix[Nts,Nts]:Rho_actor <<- Chol_to_Corr( L_Rho_actor ),
    gq> matrix[Nts,Nts]:Rho_beta <<- Chol_to_Corr( L_Rho_beta ),
    gq> matrix[Nts,Nts]:Rho_eta <<- Chol_to_Corr( L_Rho_eta )
    
  ), data = d, chains = 4, cores = 4, log_lik = TRUE )

