library(rethinking)
library(R.matlab)

fpath <- "C:\\Users\\neuro\\seadrive_root\\Emilio U\\Shared with groups\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 1\\Matlab figures\\Data\\Mice habituation.mat"
#fpath <- "C:\\Users\\Puercos\\seadrive_root\\Emilio U\\Für meine Gruppen\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 1\\Matlab figures\\Data\\Mice habituation.mat"
dat <- R.matlab::readMat(fpath)

d <- dat$mice.habituation[complete.cases(dat$mice.habituation),]
B_z <- standardize(d[,4])
d <- list(
  actor = as.integer(d[,1]),
  block_id = as.integer(d[,2]),
  tid = as.integer(2L + d[,3]*2L),
  B = B_z
)
d$tid[d$tid == 0L] = 1L

# m1.1c <- ulam(
#   alist(
#     B ~ dnorm( mu, sigma),
#     mu <- g[tid] + alpha[actor,tid] + beta[block_id,tid],
#     
#     # adaptive priors
#     vector[8]:alpha[actor] ~ multi_normal(0,Rho_actor,sigma_actor),
#     vector[8]:beta[block_id] ~ multi_normal(0,Rho_block,sigma_block),
#     
#     # fixed priors
#     g[tid] ~ dnorm(0,1),
#     sigma_actor ~ dexp(1),
#     Rho_actor ~ dlkjcorr(4),
#     sigma_block ~ dexp(1),
#     sigma ~ dexp( 1 ),
#     Rho_block ~ dlkjcorr(4)
#   ) , data=d , chains=4 , cores=4 )


m1.1nc <- ulam(
  alist(
    B ~ normal(mu, sigma),
    mu <- g[tid] + alpha[actor,tid] + beta[block_id,tid],

    #Adaptive priors
    transpars> matrix[actor,8]:alpha <-
      compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
    transpars> matrix[block_id,8]:beta <-
      compose_noncentered( sigma_beta , L_Rho_beta , z_beta ),
    matrix[8,actor]:z_actor ~ normal( 0 , 1 ),
    matrix[8,block_id]:z_beta ~ normal( 0 , 1 ),

    #Fixed priors
    g[tid] ~ normal( 0 , 1 ),
    vector[8]:sigma_actor ~ exponential( 1 ),
    cholesky_factor_corr[8]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
    vector[8]:sigma_beta ~ exponential( 1 ),
    cholesky_factor_corr[8]:L_Rho_beta ~ lkj_corr_cholesky( 2 ),
    sigma ~ exponential( 1 ),

    #Finally, compute ordinary correlation matrices from Cholesky factors
    gq> matrix[8,8]:Rho_actor <<- Chol_to_Corr( L_Rho_actor ),
    gq> matrix[8,8]:Rho_beta <<- Chol_to_Corr( L_Rho_beta ),
  ), data = d, chains = 4, cores = 4, log_lik = TRUE )

m1.2nc <- ulam(
  alist(
    B ~ normal( mu, sigma ),
    mu <- g[tid] + alpha[actor,tid] + beta[block_id,tid],
    
    # adaptive priors - non-centered
    transpars> matrix[actor,8]:alpha <-
      compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
    transpars> matrix[block_id,8]:beta <-
      compose_noncentered( sigma_block , L_Rho_block , z_block ),
    matrix[8,actor]:z_actor ~ normal( 0 , 1 ),
    matrix[8,block_id]:z_block ~ normal( 0 , 1 ),
    
    # fixed priors
    g[tid] ~ normal(0,1),
    vector[8]:sigma_actor ~ exponential(1),
    cholesky_factor_corr[8]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
    vector[8]:sigma_block ~ exponential(1),
    cholesky_factor_corr[8]:L_Rho_block ~ lkj_corr_cholesky( 2 ),
    sigma ~ exponential( 1 ),
    
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[8,8]:Rho_actor <<- Chol_to_Corr(L_Rho_actor),
    gq> matrix[8,8]:Rho_block <<- Chol_to_Corr(L_Rho_block)
  ) , data=d , chains=4 , cores=4, log_lik=TRUE )

p_post <- link( m1.2nc )
p_mu <- apply( p_post , 2 , mean )
p_ci <- apply( p_post , 2 , PI )

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