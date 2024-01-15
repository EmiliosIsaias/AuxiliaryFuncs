library(rethinking)
data(chimpanzees)
d <- chimpanzees
d$block_id <- d$block
d$treatment <- 1L + d$prosoc_left + 2L*d$condition

dat <- list(
  L = d$pulled_left,
  tid = d$treatment,
  actor = d$actor,
  block_id = as.integer(d$block_id) )

set.seed(4387510)
m14.2 <- ulam(
  alist(
    L ~ dbinom(1,p),
    logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id,tid],
    
    # adaptive priors
    vector[4]:alpha[actor] ~ multi_normal(0,Rho_actor,sigma_actor),
    vector[4]:beta[block_id] ~ multi_normal(0,Rho_block,sigma_block),
    
    # fixed priors
    g[tid] ~ dnorm(0,1),
    sigma_actor ~ dexp(1),
    Rho_actor ~ dlkjcorr(4),
    sigma_block ~ dexp(1),
    Rho_block ~ dlkjcorr(4)
  ) , data=dat , chains=4 , cores=4 )

# Non-centered parameters
set.seed(4387510)
m14.3 <- ulam(
  alist(
    L ~ binomial(1,p),
    logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id,tid],
    
    # adaptive priors - non-centered
    transpars> matrix[actor,4]:alpha <-
      compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
    transpars> matrix[block_id,4]:beta <-
      compose_noncentered( sigma_block , L_Rho_block , z_block ),
    matrix[4,actor]:z_actor ~ normal( 0 , 1 ),
    matrix[4,block_id]:z_block ~ normal( 0 , 1 ),
    
    # fixed priors
    g[tid] ~ normal(0,1),
    vector[4]:sigma_actor ~ dexp(1),
    cholesky_factor_corr[4]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
    vector[4]:sigma_block ~ dexp(1),
    cholesky_factor_corr[4]:L_Rho_block ~ lkj_corr_cholesky( 2 ),
    
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[4,4]:Rho_actor <<- Chol_to_Corr(L_Rho_actor),
    gq> matrix[4,4]:Rho_block <<- Chol_to_Corr(L_Rho_block)
  ) , data=dat , chains=4 , cores=4 , log_lik=TRUE )


neff_nc <- precis(m14.3,3,pars=c("alpha","beta"))$n_eff
neff_c <- precis(m14.2,3,pars=c("alpha","beta"))$n_eff
plot( neff_c , neff_nc , xlab="centered (default)" ,
      ylab="non-centered (cholesky)" , lwd=1.5 )
abline(a=0,b=1,lty=2)


# compute mean for each actor in each treatment
pl <- by( d$pulled_left , list( d$actor , d$treatment ) , mean )

# generate posterior predictions using link
datp <- list(
  actor=rep(1:7,each=4) ,
  tid=rep(1:4,times=7) ,
  block_id=rep(5,times=4*7) )
p_post <- link( m14.3 , data=datp )
p_mu <- apply( p_post , 2 , mean )
p_ci <- apply( p_post , 2 , PI )

# set up plot
plot( NULL , xlim=c(1,28) , ylim=c(0,1) , xlab="" ,
      ylab="proportion left lever" , xaxt="n" , yaxt="n" )
axis( 2 , at=c(0,0.5,1) , labels=c(0,0.5,1) )
abline( h=0.5 , lty=2 )
for ( j in 1:7 ) abline( v=(j-1)*4+4.5 , lwd=0.5 )
for ( j in 1:7 ) text( (j-1)*4+2.5 , 1.1 , concat("actor ",j) , xpd=TRUE )

xo <- 0.1 # offset distance to stagger raw data and predictions
# raw data
for ( j in (1:7)[-2] ) {
  lines( (j-1)*4+c(1,3)-xo , pl[j,c(1,3)] , lwd=2 , col=rangi2 )
  lines( (j-1)*4+c(2,4)-xo , pl[j,c(2,4)] , lwd=2 , col=rangi2 )
}
points( 1:28-xo , t(pl) , pch=16 , col="white" , cex=1.7 )
points( 1:28-xo , t(pl) , pch=c(1,1,16,16) , col=rangi2 , lwd=2 )

yoff <- 0.175
text( 1-xo , pl[1,1]-yoff , "R/N" , pos=1 , cex=0.8 )
text( 2-xo , pl[1,2]+yoff , "L/N" , pos=3 , cex=0.8 )
text( 3-xo , pl[1,3]-yoff , "R/P" , pos=1 , cex=0.8 )
text( 4-xo , pl[1,4]+yoff , "L/P" , pos=3 , cex=0.8 )

# posterior predictions
for ( j in (1:7)[-2] ) {
  lines( (j-1)*4+c(1,3)+xo , p_mu[(j-1)*4+c(1,3)] , lwd=2 )
  lines( (j-1)*4+c(2,4)+xo , p_mu[(j-1)*4+c(2,4)] , lwd=2 )
}
for ( i in 1:28 ) lines( c(i,i)+xo , p_ci[,i] , lwd=1 )
points( 1:28+xo , p_mu , pch=16 , col="white" , cex=1.3 )
points( 1:28+xo , p_mu , pch=c(1,1,16,16) )

# Instrumental variables
set.seed(73)
N <- 500
U_sim <- rnorm( N )
Q_sim <- sample( 1:4 , size=N , replace=TRUE )
E_sim <- rnorm( N , U_sim + Q_sim )
W_sim <- rnorm( N , U_sim + 0*E_sim )
dat_sim <- list(
  W=standardize(W_sim) ,
  E=standardize(E_sim) ,
  Q=standardize(Q_sim) )

m14.4 <- ulam(
  alist(
    W ~ dnorm( mu , sigma ),
    mu <- aW + bEW*E,
    aW ~ dnorm( 0 , 0.2 ),
    bEW ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ) , data=dat_sim , chains=4 , cores=4 )
precis( m14.4 )

m14.5 <- ulam(
  alist(
    W ~ dnorm( mu , sigma ),
    mu <- aW + bEW*E + bQW*Q,
    aW ~ dnorm( 0 , 0.2 ),
    bEW ~ dnorm( 0 , 0.5 ),
    bQW ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ) , data=dat_sim , chains=4 , cores=4 )
precis( m14.5 )

m14.6 <- ulam(
  alist(
    c(W,E) ~ multi_normal( c(muW,muE) , Rho , Sigma ),
    muW <- aW + bEW*E,
    muE <- aE + bQE*Q,
    c(aW,aE) ~ normal( 0 , 0.2 ),
    c(bEW,bQE) ~ normal( 0 , 0.5 ),
    Rho ~ lkj_corr( 2 ),
    Sigma ~ exponential( 1 )
  ), data=dat_sim , chains=4 , cores=4 )
precis( m14.6 , depth=3 )

m14.4x <- ulam( m14.4 , data=dat_sim , chains=4 , cores=4 )
m14.6x <- ulam( m14.6 , data=dat_sim , chains=4 , cores=4 )


data(KosterLeckie)

## R code 14.31
kl_data <- list(
  N = nrow(kl_dyads),
  N_households = max(kl_dyads$hidB),
  did = kl_dyads$did,
  hidA = kl_dyads$hidA,
  hidB = kl_dyads$hidB,
  giftsAB = kl_dyads$giftsAB,
  giftsBA = kl_dyads$giftsBA
)

m14.7 <- ulam(
  alist(
    giftsAB ~ poisson( lambdaAB ),
    giftsBA ~ poisson( lambdaBA ),
    log(lambdaAB) <- a + gr[hidA,1] + gr[hidB,2] + d[did,1] ,
    log(lambdaBA) <- a + gr[hidB,1] + gr[hidA,2] + d[did,2] ,
    a ~ normal(0,1),
    
    ## gr matrix of varying effects
    vector[2]:gr[N_households] ~ multi_normal(0,Rho_gr,sigma_gr),
    Rho_gr ~ lkj_corr(4),
    sigma_gr ~ exponential(1),
    
    ## dyad effects
    transpars> matrix[N,2]:d <-
      compose_noncentered( rep_vector(sigma_d,2) , L_Rho_d , z ),
    matrix[2,N]:z ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_d ~ lkj_corr_cholesky( 8 ),
    sigma_d ~ exponential(1),
    
    ## compute correlation matrix for dyads
    gq> matrix[2,2]:Rho_d <<- Chol_to_Corr( L_Rho_d )
  ), data=kl_data , chains=4 , cores=4 , iter=2000 )

post <- extract.samples( m14.7 )
g <- sapply( 1:25 , function(i) post$a + post$gr[,i,1] )
r <- sapply( 1:25 , function(i) post$a + post$gr[,i,2] )
Eg_mu <- apply( exp(g) , 2 , mean )
Er_mu <- apply( exp(r) , 2 , mean )

plot( NULL , xlim=c(0,8.6) , ylim=c(0,8.6) , xlab="generalized giving" ,
      ylab="generalized receiving" , lwd=1.5 )
abline(a=0,b=1,lty=2)

# ellipses
library(ellipse)
for ( i in 1:25 ) {
  Sigma <- cov( cbind( g[,i] , r[,i] ) )
  Mu <- c( mean(g[,i]) , mean(r[,i]) )
  for ( l in c(0.5) ) {
    el <- ellipse( Sigma , centre=Mu , level=l )
    lines( exp(el) , col=col.alpha("black",0.5) )
  }
}
# household means
points( Eg_mu , Er_mu , pch=21 , bg="white" , lwd=1.5 )