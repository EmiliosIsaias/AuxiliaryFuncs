library(rethinking)
library(R.matlab)

# fpath <- file.path("C:", "Users", "jefe_", "seadrive_root", "Emilio U", 
#                    "Für meine Gruppen", "GDrive GrohLab", "Projects", 
#                    "00 Salience", "Bayes Model Data", "Lick_data_4_R.mat", 
#                    fsep = .Platform$file.sep )
# fpath <- file.path("c:", "users", "jefe_", "seadrive_root", "emilio u",
#                    ""
#                    fsep = "\\" )

data_path <- file.path("c:", "users", "neuro", "seadrive_root", "emilio u",
                   "Shared with groups", "gdrive grohlab", "projects",
                   "00 SC", "SC Behaviour", "Figures", "Figure 2", 
                   "Matlab figures", "Data" )
fpath <- file.path( data_path, "MC_opsines_4_r.mat" )
dat <- R.matlab::readMat(fpath)

bi_mu <- dat$bi.centre[1]
bi_sig <- dat$bi.scale[1]

z_bi <- ( dat$bi[,1] - bi_mu ) / bi_sig

cc_flag = complete.cases(dat$bi)

d <- list(
  bi = z_bi[cc_flag],
  #l_lat = z_lick,
  mouse_id = as.integer( dat$mouse.id[cc_flag] ),
  tid = as.integer( dat$tid[cc_flag] ),
  Nc = as.integer( dat$Nc )
)

mMCBI <- ulam( 
  alist(
    bi ~ normal( mu , sigma ),
    mu <- g[tid] + alpha[mouse_id,tid],
    
    #Adaptive priors
    transpars> matrix[mouse_id,Nc]:alpha <-
      compose_noncentered( sigma_mouse , L_Rho_mouse , z_mouse ),
    
    matrix[Nc,mouse_id]:z_mouse ~ normal( 0 , 1 ),
    

    #Fixed priors
    g[tid] ~ normal( 0 , 1 ),
    vector[Nc]:sigma_mouse ~ exponential( 1 ),
    cholesky_factor_corr[Nc]:L_Rho_mouse ~ lkj_corr_cholesky( 2 ),
    
    sigma ~ exponential( 1 ),

    #Finally, compute ordinary correlation matrices from Cholesky factors
    gq> matrix[Nc,Nc]:Rho_mouse <<- Chol_to_Corr( L_Rho_mouse )
    
  ), data = d, chains = 4, cores = 4, log_lik = TRUE )

params <- extract.samples( mMCBI )
post <- link( mMCBI )

fpath_out <- file.path(data_path, "Opsines_Bayes.mat")
R.matlab::writeMat(fpath_out, params = params, post = post)
