library(rethinking)
library(R.matlab)

#fpath <- "C:\\Users\\neuro\\seadrive_root\\Emilio U\\Shared with groups\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 1\\Matlab figures\\Data\\Mice habituation.mat"
#fpath <- "C:\\Users\\Puercos\\seadrive_root\\Emilio U\\Für meine Gruppen\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 1\\Matlab figures\\Data\\Mice habituation.mat"
#fpath <- "C:\\Users\\jefe_\\seadrive_root\\Emilio U\\Für meine Gruppen\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 1\\Matlab figures\\Data\\Mice habituation.mat"
fpath <- "C:\\Users\\jefe_\\seadrive_root\\Emilio U\\Für meine Gruppen\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 3\\Matlab figures\\Data\\iRNs_4_stan.mat"
dat <- R.matlab::readMat(fpath)

bi_mu = mean( dat$bi, na.rm = TRUE)
bi_sd = sd( dat$bi, na.rm = TRUE)


d <- list(
  y = (dat$bi[complete.cases(dat$bi),] - bi_mu) / bi_sd,
  actor_id = as.integer( dat$mouse.id[complete.cases(dat$bi),] ),
  block_id = as.integer( dat$block.id[complete.cases(dat$bi),] ),
  tid = as.integer( dat$tid[complete.cases(dat$bi),] ),
  N = sum( complete.cases( dat$bi ) ),
  Nm = as.integer( dat$Nm[1] ),
  Nc = as.integer( dat$Nc[1] ),
  Ns = as.integer( dat$Ns[1] )
)

mIRN_MC <- ulam( alist(
    y ~ normal( mu , sigma ),
    mu <- g[tid] + alpha[actor_id,tid] + beta[block_id,tid],
    
    #Adaptive priors
    transpars> matrix[actor_id,Nc]:alpha <-
      compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
    transpars> matrix[block_id,Nc]:beta <-
      compose_noncentered( sigma_beta , L_Rho_beta , z_beta ),
    matrix[Nc,actor_id]:z_actor ~ normal( 0 , 1 ),
    matrix[Nc,block_id]:z_beta ~ normal( 0 , 1 ),

    #Fixed priors
    g[tid] ~ normal( 0 , 1 ),
    vector[Nc]:sigma_actor ~ exponential( 1 ),
    cholesky_factor_corr[Nc]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
    vector[Nc]:sigma_beta ~ exponential( 1 ),
    cholesky_factor_corr[Nc]:L_Rho_beta ~ lkj_corr_cholesky( 2 ),
    sigma ~ exponential( 1 ),

    #Finally, compute ordinary correlation matrices from Cholesky factors
    gq> matrix[Nc,Nc]:Rho_actor <<- Chol_to_Corr( L_Rho_actor ),
    gq> matrix[Nc,Nc]:Rho_beta <<- Chol_to_Corr( L_Rho_beta )
  ), data = d, chains = 4, cores = 4, log_lik = TRUE)

post <- link( mIRN_MC )
params <- extract.samples( mIRN_MC )

R.matlab::writeMat("C:/Users/jefe_/seadrive_root/Emilio U/Für meine Gruppen/GDrive GrohLab/Projects/00 SC/SC Behaviour/Figures/Figure 3/Matlab figures/Data/Bayes_iRN.mat", post = post, params = params)

xax <- seq(-0.5, 3, 0.5)
plot( NULL , xlim=c(-0.5,3) , ylim=c(0,1) , xlab="Puff intensities" ,
      ylab="Z-BI" , xaxt="n" , yaxt="n" )
abline( h=0 , lty=2 )
points(d$tid/2 - 1, p_mu, col=rangi2)


plot( NULL , xlim=c(1,28) , ylim=c(0,1) , xlab="" ,
      ylab="proportion left lever" , xaxt="n" , yaxt="n" )
axis( 2 , at=c(0,0.5,1) , labels=c(0,0.5,1) )
abline( h=0.5 , lty=2 )
for ( j in 1:7 ) abline( v=(j-1)*4+4.5 , lwd=0.5 )
for ( j in 1:7 ) text( (j-1)*4+2.5 , 1.1 , concat("actor ",j) , xpd=TRUE )