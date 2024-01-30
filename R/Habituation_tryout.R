library(rethinking)
library(R.matlab)

fpath <- "C:\\Users\\neuro\\seadrive_root\\Emilio U\\Shared with groups\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 1\\Matlab figures\\Data\\Mice habituation.mat"
dat <- R.matlab::readMat(fpath)

d <- dat$mice.habituation[complete.cases(dat$mice.habituation),]
d <- list(
  mouse = as.integer(d[,1]),
  sess = as.integer(d[,2]),
  tid = as.integer(2L + d[,3]*2L),
  B = standardize(d[,4])
)
d$tid[d$tid == 0L] = 1L

m1.2 <- ulam(
  alist(
    B ~ normal(mu, sigma),
    mu <- gama[tid] + mice[mouse,tid] + session[sess, tid],
    
    #Adaptive priors
    transpars> matrix[mouse,8]:mice <- 
      compose_noncentered(sigma_m, L_Rho_m, z_m),
    transpars> matrix[sess,8]:session <- 
      compose_noncentered(sigma_s, L_Rho_s, z_s),
    matrix[8,mouse]:z_m ~ dnorm(0, 1),
    matrix[8,sess]:z_s ~ dnorm(0, 1),
    
    #Fixed priors
    gama[tid] ~ dnorm(0, 1),
    vector[8]:sigma_m ~ dexp(1),
    vector[8]:sigma_s ~ dexp(1),
    sigma ~ dexp(1),
    cholesky_factor_corr[8]:L_Rho_m ~ lkj_corr_cholesky(2),
    cholesky_factor_corr[8]:L_Rho_s ~ lkj_corr_cholesky(2),
    
    #Finally, compute ordinary correlation matrices from Cholesky factors
    gq> matrix[8,8]:Rho_m <<- Chol_to_Corr(L_Rho_m),
    gq> matrix[8,8]:Rho_s <<- Chol_to_Corr(L_Rho_s),
  ), data = d, chains = 4, cores = 4, log_lik = TRUE )
