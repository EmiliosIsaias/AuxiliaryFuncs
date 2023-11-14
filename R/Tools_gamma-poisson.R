library(rethinking)  
data(Kline)  
d <- Kline  
d$P <- standardize( log(d$population))  
d$contact_id <- ifelse( d$contact=="high", 2L, 1L)

dat2 <- list(  
  T = d$total_tools,
  P = d$population,  
  cid = d$contact_id ) 

m12.2 <- ulam(  
  alist(  
    T ~ dgampois( lambda, phi),  
    lambda <- exp(a[cid])*P^b[cid]/ g,
    a[cid] ~ dnorm(1,1),
    b[cid] ~ dexp(1),
    g ~ dexp(1),  
    phi ~ dexp(1)  
    ), data=dat2, chains=4, log_lik=TRUE, cores = 4
  )