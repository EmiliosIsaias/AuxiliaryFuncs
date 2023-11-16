library(rethinking)
library(R.matlab)

mdata <- R.matlab::readMat( 
  "C:\\Users\\neuro\\seadrive_root\\Emilio U\\My Libraries\\My Library\\VPM-POm spatial correlation\\vpm_pom_cluster_data_for_Emilio_Oct_23.mat"
  )

td2 <- list( vpm = mdata$is.vpm[,1], 
            PF =  matrix( standardize( mdata$Puff.PSTH ) , nrow = 28 ) ,
            TC =  matrix( standardize( mdata$Touch.PSTH ) , nrow = 28 )  
            )

mth.binom <- quap(  
  alist(    
    vpm ~ dbinom(1, p),    
    logit(p) <- a,    
    a ~ dnorm(0,3)  
    ), data = td2 )

mth.p <- quap(
  alist(
    vpm ~ dbinom(1, p),    
    logit(p) <- a + sapply( 1:28, function(i) sum(PF[i,] * b) ),
    a ~ dnorm(0, 3),
    b ~ dnorm(0, 5)
  ), data = td2 )

mth.p <- quap(
  alist(
    vpm ~ dbinom(1, p),    
    logit(p) <- a + sapply( 1:100, function(i) sum(PF[i,] * b) ),
    a ~ dnorm(0, 3),
    b ~ dnorm(0, 5)
  ), data = td2 )

## R code 4.72
data(cherry_blossoms)
d <- cherry_blossoms
precis(d)

## R code 4.73
d2 <- d[ complete.cases(d$doy) , ] # complete cases on doy
library(splines)
num_knots <- 15
knot_list <- quantile( d2$year , probs=seq(0,1,length.out=num_knots) )

B <- bs(d2$year,
        knots=knot_list[-c(1,num_knots)] ,
        degree=3 , intercept=TRUE )

m4.7 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + B %*% w ,
    a ~ dnorm(100,10),
    w ~ dnorm(0,10),
    sigma ~ dexp(1)
  ), data=list( D=d2$doy , B=B ) , 
  start=list( w=rnorm(ncol(B), 0, 10) ) )

post <- extract.samples( m4.7 )
w <- apply( post$w , 2 , mean )
plot( NULL , xlim=range(d2$year) , ylim=c(-6,6) ,
      xlab="year" , ylab="basis * weight" )
for ( i in 1:ncol(B) ) lines( d2$year , w[i]*B[,i] )