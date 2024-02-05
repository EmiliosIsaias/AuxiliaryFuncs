library(rethinking)
library(R.matlab)

fpath <- "C:\\Users\\jefe_\\seadrive_root\\Emilio U\\Für meine Gruppen\\GDrive GrohLab\\Projects\\00 Jitter\\Bayes_psth_test_1_neu_spkTms.mat"

dat <- R.matlab::readMat(fpath)

d <- list(
  y = dat$unit.spike.times[,3] ,
  trig_id = as.integer( dat$unit.spike.times[,2] ),
  cond_id = as.integer( dat$unit.spike.times[,1] )
)

m1.nc <- ulam(
  alist(
    y ~ dzipois( p , lambda ),
    mu <- g[cond_id] + beta[trig_id,cond_id],
    
    # adaptive priors - non-centered
    transpars> matrix[trig_id,6]:beta <-
      compose_noncentered( sigma_block , L_Rho_block , z_block ),
    matrix[6,trig_id]:z_block ~ normal( 0 , 1 ),
    
    # fixed priors
    g[cond_id] ~ normal(0,1),
    vector[6]:sigma_block ~ dexp(1),
    cholesky_factor_corr[6]:L_Rho_block ~ lkj_corr_cholesky( 2 ),
    sigma ~ dexp( 1 ),
    
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[6,6]:Rho_block <<- Chol_to_Corr(L_Rho_block)
  ) , data=d , chains=4 , cores=4 , log_lik=TRUE )

p_post <- link(m1.nc, 
               data = list(
                 trig_id = rep(1:40, times = 6),
                 cond_id = rep(1:6, each = 40)
               ) )
post <- extract.samples(m1.nc)

R.matlab::writeMat("C:\\Users\\jefe_\\seadrive_root\\Emilio U\\Für meine Gruppen\\GDrive GrohLab\\Projects\\00 Jitter\\Post_1neu_chimp.mat", post = post, p_post = p_post )
