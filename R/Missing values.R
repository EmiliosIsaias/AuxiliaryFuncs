## R code 15.2
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce

# points
plot( d$Divorce ~ d$MedianAgeMarriage , ylim=c(4,15) ,
      xlab="Median age marriage" , ylab="Divorce rate" )

# standard errors
for ( i in 1:nrow(d) ) {
  ci <- d$Divorce[i] + c(-1,1)*d$Divorce.SE[i]
  x <- d$MedianAgeMarriage[i]
  lines( c(x,x) , ci )
}

## R code 15.3
dlist <- list(
  D_obs = standardize( d$Divorce ),
  D_sd = d$Divorce.SE / sd( d$Divorce ),
  M = standardize( d$Marriage ),
  A = standardize( d$MedianAgeMarriage ),
  N = nrow(d)
)


m15.1 <- ulam(
  alist(
    D_obs ~ dnorm( D_true , D_sd ),
    vector[N]:D_true ~ dnorm( mu , sigma ),
    mu <- a + bA*A + bM*M,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ) , data=dlist , chains=4 , cores=4 )

## R code 15.4
precis( m15.1 , depth=2 )

## R code 15.8
N <- 100
S <- rnorm( N )
H <- rbinom( N , size=10 , inv_logit(S) )

## R code 15.9
D <- rbern( N ) # dogs completely random
Hm <- H
Hm[D==1] <- NA

dlist <- list(
  S = S,
  H = Hm
)

m15.0 <- ulam(
  alist(
    H ~ dbinom( 10,  p ),
    logit(p) <- a + bS*S,
    a ~ dnorm(0,1),
    bS ~ dnorm(0,1)
  ), data = dlist, chains = 4  )

## R code 15.10 dog eats good student's homework
D <- ifelse( S > 0 , 1 , 0 )
Hm <- H
Hm[D==1] <- NA

## R code 15.11
set.seed(501)
N <- 1000
X <- rnorm(N)
S <- rnorm(N)
H <- rbinom( N , size=10 , inv_logit( 2 + S - 2*X ) )
D <- ifelse( X > 1 , 1 , 0 )
Hm <- H
Hm[D==1] <- NA

## R code 15.12
dat_list <- list(
  H = H,
  S = S )

m15.3 <- ulam(
  alist(
    H ~ dbinom( 10 , p ),
    logit(p) <- a + bS*S,
    a ~ normal( 0 , 1 ),
    bS ~ normal( 0 , 0.5 )
  ), data=dat_list , chains=4 )
precis( m15.3 )

## R code 15.13
dat_list0 <- list( H = H[D==0] , S = S[D==0] )

m15.4 <- ulam(
  alist(
    H ~ binomial( 10 , p ),
    logit(p) <- a + bS*S,
    a ~ normal( 0 , 1 ),
    bS ~ normal( 0 , 0.5 )
  ), data=dat_list0 , chains=4 )
precis( m15.4 )

## R code 15.14
D <- ifelse( abs(X) < 1 , 1 , 0 )

## R code 15.15
N <- 100
S <- rnorm(N)
H <- rbinom( N , size=10 , inv_logit(S) )
D <- ifelse( H < 5 , 1 , 0 )
Hm <- H; Hm[D==1] <- NA

## R code 15.16
library(rethinking)
data(milk)
d <- milk
d$neocortex.prop <- d$neocortex.perc / 100
d$logmass <- log(d$mass)
dat_list <- list(
  K = standardize( d$kcal.per.g ),
  B = standardize( d$neocortex.prop ),
  M = standardize( d$logmass ) )

## R code 15.17
m15.5 <- ulam(
  alist(
    K ~ dnorm( mu , sigma ),
    mu <- a + bB*B + bM*M,
    B ~ dnorm( nu , sigma_B ),
    c(a,nu) ~ dnorm( 0 , 0.5 ),
    c(bB,bM) ~ dnorm( 0, 0.5 ),
    sigma_B ~ dexp( 1 ),
    sigma ~ dexp( 1 )
  ) , data=dat_list , chains=4 , cores=4 )
## R code 15.18
precis( m15.5 , depth=2 )

## R code 15.19
obs_idx <- which( !is.na(d$neocortex.prop) )
dat_list_obs <- list(
  K = dat_list$K[obs_idx],
  B = dat_list$B[obs_idx],
  M = dat_list$M[obs_idx] )
m15.6 <- ulam(
  alist(
    K ~ dnorm( mu , sigma ),
    mu <- a + bB*B + bM*M,
    B ~ dnorm( nu , sigma_B ),
    c(a,nu) ~ dnorm( 0 , 0.5 ),
    c(bB,bM) ~ dnorm( 0, 0.5 ),
    sigma_B ~ dexp( 1 ),
    sigma ~ dexp( 1 )
  ) , data=dat_list_obs , chains=4 , cores=4 )
precis( m15.6 )

## R code 15.20
plot( coeftab(m15.5,m15.6) , pars=c("bB","bM") )

## R code 15.21
post <- extract.samples( m15.5 )
B_impute_mu <- apply( post$B_impute , 2 , mean )
B_impute_ci <- apply( post$B_impute , 2 , PI )

# B vs K
plot( dat_list$B , dat_list$K , pch=16 , col=rangi2 ,
      xlab="neocortex percent (std)" , ylab="kcal milk (std)" )
miss_idx <- which( is.na(dat_list$B) )
Ki <- dat_list$K[miss_idx]
points( B_impute_mu , Ki )
for ( i in 1:12 ) lines( B_impute_ci[,i] , rep(Ki[i],2) )

# M vs B
plot( dat_list$M , dat_list$B , pch=16 , col=rangi2 ,
      ylab="neocortex percent (std)" , xlab="log body mass (std)" )
Mi <- dat_list$M[miss_idx]
points( Mi , B_impute_mu )
for ( i in 1:12 ) lines( rep(Mi[i],2) , B_impute_ci[,i] )

## R code 15.22
m15.7 <- ulam(
  alist(
    # K as function of B and M
    K ~ dnorm( mu , sigma ),
    mu <- a + bB*B_merge + bM*M,
    
    # M and B correlation
    MB ~ multi_normal( c(muM,muB) , Rho_BM , Sigma_BM ),
    matrix[29,2]:MB <<- append_col( M , B_merge ),
    
    # define B_merge as mix of observed and imputed values
    vector[29]:B_merge <- merge_missing( B , B_impute ),
    
    # priors
    c(a,muB,muM) ~ dnorm( 0 , 0.5 ),
    c(bB,bM) ~ dnorm( 0, 0.5 ),
    sigma ~ dexp( 1 ),
    Rho_BM ~ lkj_corr(2),
    Sigma_BM ~ dexp(1)
  ) , data=dat_list , chains=4 , cores=4 )
precis( m15.7 , depth=3 , pars=c("bM","bB","Rho_BM" ) )

## R code 15.23
B_missidx <- which( is.na( dat_list$B ) )

## R code 15.24
data(Moralizing_gods)
str(Moralizing_gods)

## R code 15.25
table( Moralizing_gods$moralizing_gods , useNA="always" )

## R code 15.26
symbol <- ifelse( Moralizing_gods$moralizing_gods==1 , 16 , 1 )
symbol <- ifelse( is.na(Moralizing_gods$moralizing_gods) , 4 , symbol )
color <- ifelse( is.na(Moralizing_gods$moralizing_gods) , "black" , rangi2 )
plot( Moralizing_gods$year , Moralizing_gods$population , pch=symbol ,
      col=color , xlab="Time (year)" , ylab="Population size" , lwd=1.5 )

## R code 15.27
with( Moralizing_gods ,
      table( gods=moralizing_gods , literacy=writing , useNA="always" ) )

## R code 15.28
haw <- which( Moralizing_gods$polity=="Big Island Hawaii" )
columns <- c("year","writing","moralizing_gods")
t( Moralizing_gods[ haw , columns ] )

## R code 15.29
set.seed(9)
N_houses <- 100L
alpha <- 5
beta <- (-3)
k <- 0.5
r <- 0.2
cat <- rbern( N_houses , k )
notes <- rpois( N_houses , alpha + beta*cat )
R_C <- rbern( N_houses , r )
cat_obs <- cat
cat_obs[R_C==1] <- (-9L)
dat <- list(
  notes = notes,
  cat = cat_obs,
  RC = R_C,
  N = as.integer(N_houses) )

## R code 15.30
m15.8 <- ulam(
  alist(
    # singing bird model
    ## cat known present/absent:
    notes|RC==0 ~ poisson( lambda ),
    log(lambda) <- a + b*cat,
    ## cat NA:
    notes|RC==1 ~ custom( log_sum_exp(
      log(k) + poisson_lpmf( notes | exp(a + b) ),
      log(1-k) + poisson_lpmf( notes | exp(a) )
    ) ),
    
    # priors
    a ~ normal(0,1),
    b ~ normal(0,0.5),
    
    # sneaking cat model
    cat|RC==0 ~ bernoulli(k),
    k ~ beta(2,2)
  ), data=dat , chains=4 , cores=4 )

## R code 15.31
m15.9 <- ulam(
  alist(
    # singing bird model
    notes|RC==0 ~ poisson( lambda ),
    notes|RC==1 ~ custom( log_sum_exp(
      log(k) + poisson_lpmf( notes | exp(a + b) ),
      log(1-k) + poisson_lpmf( notes | exp(a) )
    ) ),
    log(lambda) <- a + b*cat,
    a ~ normal(0,1),
    b ~ normal(0,0.5),
    
    # sneaking cat model
    cat|RC==0 ~ bernoulli(k),
    k ~ beta(2,2),
    
    # imputed values
    gq> vector[N]:PrC1 <- exp(lpC1)/(exp(lpC1)+exp(lpC0)),
    gq> vector[N]:lpC1 <- log(k) + poisson_lpmf( notes[i] | exp(a+b) ),
    gq> vector[N]:lpC0 <- log(1-k) + poisson_lpmf( notes[i] | exp(a) )
  ), data=dat , chains=4 , cores=4 )

## R code 15.32
set.seed(100)
x <- c( rnorm(10) , NA )
y <- c( rnorm(10,x) , 100 )
d <- list(x=x,y=y)

## R code 15.33
## R code 15.34
## R code 15.35
## R code 15.36
## R code 15.37
## R code 15.38
## R code 15.39
## R code 16.1
library(rethinking)
data(Howell1)
d <- Howell1

# scale observed variables
d$w <- d$weight / mean(d$weight)
d$h <- d$height / mean(d$height)

## R code 16.2
m16.1 <- ulam(
  alist(
    w ~ dlnorm( mu , sigma ),
    exp(mu) <- 3.141593 * k * p^2 * h^3,
    p ~ beta( 2 , 18 ),
    k ~ exponential( 0.5 ),
    sigma ~ exponential( 1 )
  ), data=d , chains=4 , cores=4 )

## R code 16.3
h_seq <- seq( from=0 , to=max(d$h) , length.out=30 )
w_sim <- sim( m16.1 , data=list(h=h_seq) )
mu_mean <- apply( w_sim , 2 , mean )
w_CI <- apply( w_sim , 2 , PI )
plot( d$h , d$w , xlim=c(0,max(d$h)) , ylim=c(0,max(d$w)) , col=rangi2 ,
      lwd=2 , xlab="height (scaled)" , ylab="weight (scaled)" )
lines( h_seq , mu_mean )
shade( w_CI , h_seq )

## R code 16.4
library(rethinking)
data(Boxes)
precis(Boxes)

## R code 16.5
table( Boxes$y ) / length( Boxes$y )

## R code 16.6
set.seed(7)
N <- 30 # number of children

# half are random
# sample from 1,2,3 at random for each
y1 <- sample( 1:3 , size=N/2 , replace=TRUE )

# half follow majority
y2 <- rep( 2 , N/2 )

# combine and shuffle y1 and y2
y <- sample( c(y1,y2) )

# count the 2s
sum(y==2)/N

## R code 16.7
data(Boxes_model)
cat(Boxes_model)

## R code 16.8
# prep data
dat_list <- list(
  N = nrow(Boxes),
  y = Boxes$y,
  majority_first = Boxes$majority_first )

# run the sampler
Boxes_model2 <- "
data{
    int N;
    array[N] int y;
    array[N] int majority_first;
}
parameters{
    simplex[5] p;
}
model{
    vector[5] phi;
    
    // prior
    p ~ dirichlet( rep_vector(4,5) );
    
    // probability of data
    for ( i in 1:N ) {
        if ( y[i]==2 ) phi[1]=1; else phi[1]=0; // majority
        if ( y[i]==3 ) phi[2]=1; else phi[2]=0; // minority
        if ( y[i]==1 ) phi[3]=1; else phi[3]=0; // maverick
        phi[4]=1.0/3.0;                         // random
        if ( majority_first[i]==1 )             // follow first
            if ( y[i]==2 ) phi[5]=1; else phi[5]=0;
        else
            if ( y[i]==3 ) phi[5]=1; else phi[5]=0;
        
        // compute log( p_s * Pr(y_i|s )
        for ( j in 1:5 ) phi[j] = log(p[j]) + log(phi[j]);
        // compute average log-probability of y_i
        target += log_sum_exp( phi );
    }
}"
# Original Boxes_model has an outdated definition syntax.
#m16.2 <- stan( model_code=Boxes_model , data=dat_list , chains=3 , cores=3 )
m16.2 <- stan( model_code=Boxes_model2 , data=dat_list , chains=3 , cores=3 )

# show marginal posterior for p
p_labels <- c("1 Majority","2 Minority","3 Maverick","4 Random",
              "5 Follow First")
plot( precis(m16.2,2) , labels=p_labels )

## R code 16.9
library(rethinking)
data(Panda_nuts)

## R code 16.10
N <- 1e4
phi <- rlnorm( N , log(1) , 0.1 )
k <- rlnorm( N , log(2), 0.25 )
theta <- rlnorm( N , log(5) , 0.25 )

# relative grow curve
plot( NULL , xlim=c(0,1.5) , ylim=c(0,1) , xaxt="n" , xlab="age" ,
      ylab="body mass" )
at <- c(0,0.25,0.5,0.75,1,1.25,1.5)
axis( 1 , at=at , labels=round(at*max(Panda_nuts$age)) )
for ( i in 1:20 ) curve( (1-exp(-k[i]*x)) , add=TRUE , col=grau() , lwd=1.5 )

# implied rate of nut opening curve
plot( NULL , xlim=c(0,1.5) , ylim=c(0,1.2) , xaxt="n" , xlab="age" ,
      ylab="nuts per second" )
at <- c(0,0.25,0.5,0.75,1,1.25,1.5)
axis( 1 , at=at , labels=round(at*max(Panda_nuts$age)) )
for ( i in 1:20 ) curve( phi[i]*(1-exp(-k[i]*x))^theta[i] , add=TRUE ,
                         col=grau() , lwd=1.5 )

## R code 16.11
dat_list <- list(
  n = as.integer( Panda_nuts$nuts_opened ),
  age = Panda_nuts$age / max(Panda_nuts$age),
  seconds = Panda_nuts$seconds )

m16.4 <- ulam(
  alist(
    n ~ poisson( lambda ),
    lambda <- seconds*phi*(1-exp(-k*age))^theta,
    phi ~ lognormal( log(1) , 0.1 ),
    k ~ lognormal( log(2) , 0.25 ),
    theta ~ lognormal( log(5) , 0.25 )
  ), data=dat_list , chains=4 , cores = 3)

## R code 16.12
post <- extract.samples(m16.4)
plot( NULL , xlim=c(0,1) , ylim=c(0,1.5) , xlab="age" ,
      ylab="nuts per second" , xaxt="n" )
at <- c(0,0.25,0.5,0.75,1,1.25,1.5)
axis( 1 , at=at , labels=round(at*max(Panda_nuts$age)) )

# raw data
pts <- dat_list$n / dat_list$seconds
point_size <- normalize( dat_list$seconds )
points( jitter(dat_list$age) , pts , col=rangi2 , lwd=2 , cex=point_size*3 )

# 30 posterior curves
for ( i in 1:30 ) with( post , 
                        curve( phi[i]*(1-exp(-k[i]*x))^theta[i] , add=TRUE , col=grau() ) )

## R code 16.13
library(rethinking)
data(Lynx_Hare)
plot( 1:21 , Lynx_Hare[,3] , ylim=c(0,90) , xlab="year" ,
      ylab="thousands of pelts" , xaxt="n" , type="l" , lwd=1.5 )
at <- c(1,11,21)
axis( 1 , at=at , labels=Lynx_Hare$Year[at] )
lines( 1:21 , Lynx_Hare[,2] , lwd=1.5 , col=rangi2 )
points( 1:21 , Lynx_Hare[,3] , bg="black" , col="white" , pch=21 , cex=1.4 )
points( 1:21 , Lynx_Hare[,2] , bg=rangi2 , col="white" , pch=21 , cex=1.4 )
text( 17 , 80 , "Lepus" , pos=2 )
text( 19 , 50 , "Lynx" , pos=2 , col=rangi2 )

## R code 16.14
sim_lynx_hare <- function( n_steps , init , theta , dt=0.002 ) {
  L <- rep(NA,n_steps)
  H <- rep(NA,n_steps)
  L[1] <- init[1]
  H[1] <- init[2]
  for ( i in 2:n_steps ) {
    H[i] <- H[i-1] + dt*H[i-1]*( theta[1] - theta[2]*L[i-1] )
    L[i] <- L[i-1] + dt*L[i-1]*( theta[3]*H[i-1] - theta[4] )
  }
  return( cbind(L,H) )
}

## R code 16.15
theta <- c( 0.5 , 0.05 , 0.025 , 0.5 )
#z <- sim_lynx_hare( 1e4 , as.numeric(Lynx_Hare[1,2:3]) , theta )
z <- sim_lynx_hare( 1e4 , c(1,1) , theta )

plot( z[,2] , type="l" , ylim=c(0,max(z[,2])) , lwd=2 , xaxt="n" ,
      ylab="number (thousands)" , xlab="" )
lines( z[,1] , col=rangi2 , lwd=2 )
mtext( "time" , 1 )

## R code 16.16
N <- 1e4
Ht <- 1e4
p <- rbeta(N,2,18)
h <- rbinom( N , size=Ht , prob=p )
h <- round( h/1000 , 2 )
dens( h , xlab="thousand of pelts" , lwd=2 )

## R code 16.17
data(Lynx_Hare_model)
cat(Lynx_Hare_model)

## R code 16.18
dat_list <- list(
  N = nrow(Lynx_Hare),
  pelts = Lynx_Hare[,2:3] )

Lynx_Hare_model2 <- "
functions {
  array[] real dpop_dt( real t,                 // time
                array[] real pop_init,          // initial state {lynx, hares}
                array[] real theta,             // parameters
                array[] real x_r, array[] int x_i) {  // unused
    real L = pop_init[1];
    real H = pop_init[2];
    real bh = theta[1];
    real mh = theta[2];
    real ml = theta[3];
    real bl = theta[4];
    // differential equations
    real dH_dt = (bh - mh * L) * H;
    real dL_dt = (bl * H - ml) * L;
    return { dL_dt , dH_dt };
  }
}
data {
  int<lower=0> N;              // number of measurement times
  array[N, 2] real<lower=0> pelts;    // measured populations
}
transformed data{
  array [N-1] real times_measured;    // N-1 because first time is initial state
  for ( i in 2:N ) times_measured[i-1] = i;
}
parameters {
  array[4] real<lower=0> theta;      // { bh, mh, ml, bl }
  array[2] real<lower=0> pop_init;   // initial population state
  array[2] real<lower=0> sigma;      // measurement errors
  array[2] real<lower=0,upper=1> p;  // trap rate
}
transformed parameters {
  array[N, 2] real pop;
  pop[1,1] = pop_init[1];
  pop[1,2] = pop_init[2];
  pop[2:N,1:2] = integrate_ode_rk45(
    dpop_dt, pop_init, 0, times_measured, theta,
    rep_array(0.0, 0), rep_array(0, 0),
    1e-5, 1e-3, 5e2);
}
model {
  // priors
  theta[{1,3}] ~ normal( 1 , 0.5 );    // bh,ml
  theta[{2,4}] ~ normal( 0.05, 0.05 ); // mh,bl
  sigma ~ exponential( 1 );
  pop_init ~ lognormal( log(10) , 1 );
  p ~ beta(40,200);
  // observation model
  // connect latent population state to observed pelts
  for ( t in 1:N )
    for ( k in 1:2 )
      pelts[t,k] ~ lognormal( log(pop[t,k]*p[k]) , sigma[k] );
}
generated quantities {
  array[N, 2] real pelts_pred;
  for ( t in 1:N )
    for ( k in 1:2 )
      pelts_pred[t,k] = lognormal_rng( log(pop[t,k]*p[k]) , sigma[k] );
}"
# Same case with this STAN code: outdated.
#m16.5 <- stan( model_code=Lynx_Hare_model , data=dat_list , chains=3 ,
#               cores=3 , control=list( adapt_delta=0.95 ) )
m16.5 <- stan( model_code=Lynx_Hare_model2 , data=dat_list , chains=3 ,
               cores=3 , control=list( adapt_delta=0.95 ) )

## R code 16.19
post <- extract.samples(m16.5)
pelts <- dat_list$pelts
plot( 1:21 , pelts[,2] , pch=16 , ylim=c(0,120) , xlab="year" ,
      ylab="thousands of pelts" , xaxt="n" )
at <- c(1,11,21)
axis( 1 , at=at , labels=Lynx_Hare$Year[at] )
points( 1:21 , pelts[,1] , col=rangi2 , pch=16 )
# 21 time series from posterior
for ( s in 1:21 ) {
  lines( 1:21 , post$pelts_pred[s,,2] , col=col.alpha("black",0.2) , lwd=2 )
  lines( 1:21 , post$pelts_pred[s,,1] , col=col.alpha(rangi2,0.3) , lwd=2 )
}
# text labels
text( 17 , 90 , "Lepus" , pos=2 )
text( 19 , 50 , "Lynx" , pos=2 , col=rangi2 )

## R code 16.20
plot( NULL , pch=16 , xlim=c(1,21) , ylim=c(0,500) , xlab="year" ,
      ylab="thousands of animals" , xaxt="n" )
at <- c(1,11,21)
axis( 1 , at=at , labels=Lynx_Hare$Year[at] )
for ( s in 1:21 ) {
  lines( 1:21 , post$pop[s,,2] , col=col.alpha("black",0.2) , lwd=2 )
  lines( 1:21 , post$pop[s,,1] , col=col.alpha(rangi2,0.4) , lwd=2 )
}

## R code 16.21
data(Lynx_Hare)
dat_ar1 <- list(
  L = Lynx_Hare$Lynx[2:21],
  L_lag1 = Lynx_Hare$Lynx[1:20],
  H = Lynx_Hare$Hare[2:21],
  H_lag1 = Lynx_Hare$Hare[1:20] )