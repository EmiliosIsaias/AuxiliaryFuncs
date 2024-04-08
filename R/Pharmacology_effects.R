library(rethinking)
library(R.matlab)

fpath <- "C:\\Users\\neuro\\seadrive_root\\Emilio U\\Shared with groups\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 1\\Matlab figures\\Data\\PTX_MUS_4_Stan.mat"

dat <- R.matlab::readMat(fpath)

bi_centre <- dat$bi.centre[1,]
bi_scale <- dat$bi.scale[1,]

Nt <- max( dat$tid )

d <- list(
  bi = ( dat$bi[,1] - bi_centre ) / bi_scale,
  tid = as.integer( dat$tid ),
  Nt = as.integer( Nt ),
  mouse_id = as.integer( dat$mouse.id[,1] )
)

mPharma <- ulam(
  alist(
    bi ~ normal( mu, sigma ),
    mu <- g[tid] + alpha[mouse_id, tid],
    
    transpars> matrix[mouse_id,Nt]:alpha <-
      compose_noncentered( sigma_mouse , L_Rho_mouse , z_mouse ),
    matrix[Nt,mouse_id]:z_mouse ~ normal( 0 , 1 ),
    
    # fixed priors
    g[tid] ~ normal(0,1),
    vector[Nt]:sigma_mouse ~ exponential(1),
    cholesky_factor_corr[Nt]:L_Rho_mouse ~ lkj_corr_cholesky( 2 ),
    
    sigma ~ exponential( 1 ),
    
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[Nt,Nt]:Rho_mouse <<- Chol_to_Corr(L_Rho_mouse)
  ) , data = d , chains = 4 , cores = 4 , log_lik = TRUE )

post <- link( mPharma )
params <- extract.samples( mPharma )

R.matlab::writeMat("C:\\Users\\neuro\\seadrive_root\\Emilio U\\Shared with groups\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 1\\Matlab figures\\Data\\PTX and muscimol Bayes.mat", 
                   post = post, params = params, bi_centre = bi_centre, bi_scale = bi_scale)
