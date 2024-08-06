library(rethinking)
library(R.matlab)

#fpath <- "C:\\Users\\neuro\\seadrive_root\\Emilio U\\Shared with groups\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 1\\Matlab figures\\Data\\Mice habituation.mat"
#fpath <- "C:\\Users\\Puercos\\seadrive_root\\Emilio U\\Für meine Gruppen\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 1\\Matlab figures\\Data\\Mice habituation.mat"
#fpath <- "C:\\Users\\jefe_\\seadrive_root\\Emilio U\\Für meine Gruppen\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 1\\Matlab figures\\Data\\Mice habituation.mat"
#dat <- R.matlab::readMat(fpath)

d <- dat$mice.habituation[complete.cases(dat$mice.habituation),]
B_z <- standardize(d[,4])
d <- list(
  actor = as.integer(d[,1]),
  block_id = as.integer(d[,2]),
  tid = as.integer(2L + d[,3]*2L),
  B = B_z
)

dat <- R.matlab::readMat("C:/Users/neuro/seadrive_root/Emilio U/Shared with groups/GDrive GrohLab/Projects/00 SC/SC Behaviour/Figures/Figure 1/Matlab figures/Data/Habituation_data.mat")
dat <- dat$hab.table
d <- list(
  mouse = as.integer(dat[,1]),
  day = as.integer(dat[,2]),
  intensity = dat[dat[,1]>3,3],
  B = standardize(dat[dat[,1]>3,4]),
  N = sum(dat[,1]>3)
  )
centre <- mean(dat[dat[,1]>3,4])
scale <- sd(dat[dat[,1]>3,4])

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

m1lr <- ulam( alist( 
  B ~ dnorm(mu, sigma),
  mu <- alpha + beta * intensity,
  
  beta ~ dnorm(0.4, 0.7),
  alpha ~ dnorm(0, 1),
  sigma ~ dexp(1)
  
  ), data=d, chains=4, cores=4 )

post <- extract.samples(m1lr)

alpha <- mean(post$alpha)
beta <- mean(post$beta)
sigma <- mean(post$sigma)

puff.seq <- seq( from=-0.5 , to=3.5 , by=0.2 )
mu.link <- function(x) (post$alpha + post$beta * x)*scale + centre
mu <- sapply( puff.seq, mu.link )
mu.PI <- apply(mu, 2, PI, 0.95)


plot( centre + scale*d$B ~ d$intensity )
curve( (alpha + beta * x)*scale + centre, add=TRUE)
shade(mu.PI, puff.seq, col = col.alpha("black", alpha = 1/3))



m1.1nc <- ulam( alist(
    B ~ normal( mu , sigma ),
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
    gq> matrix[8,8]:Rho_beta <<- Chol_to_Corr( L_Rho_beta )
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

m1.3nc <- ulam(
  alist(
    B ~ normal( mu, sigma ),
    mu <- alpha[mouse,day] + beta[mouse,day] * intensity,
    
    # adaptive priors - non-centered
    transpars> matrix[mouse,10]:alpha <-
      compose_noncentered( sigma_intercept , L_Rho_intercept , z_inter ),
    transpars> matrix[day,10]:beta <-
      compose_noncentered( sigma_slope , L_Rho_slope , z_slope ),
    matrix[10,mouse]:z_inter ~ normal( 0 , 1 ),
    matrix[10,day]:z_slope ~ normal( 0 , 1 ),
    
    # fixed priors
    vector[10]:sigma_intercept ~ exponential(1),
    cholesky_factor_corr[10]:L_Rho_intercept ~ lkj_corr_cholesky( 2 ),
    vector[10]:sigma_slope ~ exponential(1),
    cholesky_factor_corr[10]:L_Rho_slope ~ lkj_corr_cholesky( 2 ),
    sigma ~ exponential( 1 ),
    
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[10,10]:Rho_intercept <<- Chol_to_Corr(L_Rho_intercept),
    gq> matrix[10,10]:Rho_slope <<- Chol_to_Corr(L_Rho_slope)
  ) , data=d , chains=4 , cores=4, log_lik=TRUE )

p_post <- link( m1.3nc )
p_mu <- apply( p_post , 2 , mean )
p_ci <- apply( p_post , 2 , PI )

xax <- seq(-0.5, 3, 0.5)
plot( NULL , xlim=c(-0.5,3) , ylim=c(-1,1) , xlab="Puff intensities" ,
      ylab="Z-BI" , xaxt="n" , yaxt="n" )
abline( h=0 , lty=2 )
points(d$intensity, p_mu, col=rangi2)


plot( NULL , xlim=c(1,28) , ylim=c(0,1) , xlab="" ,
      ylab="proportion left lever" , xaxt="n" , yaxt="n" )
axis( 2 , at=c(0,0.5,1) , labels=c(0,0.5,1) )
abline( h=0.5 , lty=2 )
for ( j in 1:7 ) abline( v=(j-1)*4+4.5 , lwd=0.5 )
for ( j in 1:7 ) text( (j-1)*4+2.5 , 1.1 , concat("mouse ",j) , xpd=TRUE )