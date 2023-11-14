library(rethinking)  
data(UCBadmit)  
d <- UCBadmit
d$gid <- ifelse( d$applicant.gender=="male", 1L, 2L)
dat <- list( A=d$admit, N=d$applications, gid=d$gid)
m12.1 <- ulam(  
  alist(  
    A ~ dbetabinom( N, pbar, theta),  
    logit(pbar) <- a[gid],
    a[gid] ~ dnorm( 0, 1.5),  
    transpars> theta <<- phi + 2.0,  
    phi ~ dexp(1)  
    ), 
  data=dat, chains=4, cores = 4
  ) 

post <- extract.samples(m12.1)
post$da <- post$a[,1] - post$a[,2]
precis(post, depth = 2)

gid <- 2
#draw posterior mean beta distribution  
curve( dbeta2(x,mean(logistic(post$a[,gid])),mean(post$theta)), from=0, to=1,  
       ylab="Density", xlab="probability admit", ylim=c(0,3), lwd=2)

#draw 50 beta distributions sampled from posterior 
for (i in 1:50) {  
  p <- logistic( post$a[i,gid])
  theta <- post$theta[i]  
  curve( dbeta2(x,p,theta), add=TRUE, col=col.alpha("black",0.2))  
  }  
mtext( "distribution of female admission rates") 