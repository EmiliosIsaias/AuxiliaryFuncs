library(rethinking)
library(R.matlab)

#mdata <- R.matlab::readMat( "C:\\Users\\neuro\\seadrive_root\\Emilio U\\My Libraries\\My Library\\VPM-POm spatial correlation\\vpm_pom_cluster_data_for_Emilio_Oct_23.mat" )
mdata <- R.matlab::readMat( 
  "C:\\Users\\jefe_\\seadrive_root\\Emilio U\\Meine Bibliotheken\\My Library\\VPM-POm spatial correlation\\vpm_pom_cluster_data_for_Emilio_Oct_23.mat"
)

#mdata2 <- R.matlab::readMat(
#  "C:\\Users\\neuro\\seadrive_root\\Emilio U\\My Libraries\\My Library\\VPM-POm spatial correlation\\vpm_pom_cluster_data_for_Emilio.mat"
#)

# td <- list( vpm = mdata$is.vpm[,1],
#             PF =  matrix( standardize( mdata$Puff.PSTH ) , nrow = nrow(mdata$Puff.PSTH) ) ,
#             TC =  matrix( standardize( mdata$Touch.PSTH ) , nrow = nrow(mdata$Touch.PSTH) )  )

td <- list( vpm = mdata$is.vpm[,1],
            PF =  mdata$Puff.PSTH / sapply(1:nrow(mdata$Puff.PSTH), 
                                           function(i) max( mdata$Puff.PSTH[i,] ) ) ,
            TC =  mdata$Touch.PSTH / sapply(1:nrow(mdata$Touch.PSTH), 
                                            function(i) max( mdata$Touch.PSTH[i,] ) ) )


#td2 <- list( vpm = mdata2$is.vpm[,1],
#             PF =  matrix( standardize( mdata2$Puff.PSTH ) , nrow = nrow(mdata2$Puff.PSTH) ) ,
#             TC =  matrix( standardize( mdata2$Touch.PSTH ) , nrow = nrow(mdata2$Touch.PSTH) )  )

mth.binom <- quap(
  alist(
    vpm ~ dbinom(1, p),
    logit(p) <- a,
    a ~ dnorm(0,3)
    ), data = td2 )

mth.p <- quap(
  alist(
    vpm ~ dbinom(1, p),
    logit(p) <- A + PF %*% B,
    A ~ dnorm(0, 3),
    B ~ dnorm(0, 5) ),
  data = td2,   
  start = list( B = rnorm(ncol(td2$PF), 0, 5) ) )

mth.t <- quap(
  alist(
    vpm ~ dbinom(1, p),
    logit(p) <- A + TC %*% w,
    A ~ dnorm(0.57, 0.3),
    w ~ dnorm(0, 15) ), 
  data = td2,
  start = list( w = rep(0, ncol(td2$PF) ) ) )

mth.pt <- quap(
  alist(
    vpm ~ dbinom(1, p),
    logit(p) <- A + PF %*% B + TC %*% W,
    A ~ dnorm(0, 3),
    B ~ dnorm(0, 5),
    W ~ dnorm(0, 5) ), 
  data = td,   
  start = list( B = rep(0, ncol(td$PF) ),
                W = rep(0, ncol(td$TC) ) ) )

mth.pta <- quap(
  alist(
    vpm ~ dbinom(1, p),
    logit(p) <- PF %*% B + TC %*% W,
    B ~ dnorm(0, 5),
    W ~ dnorm(0, 5) ), 
  data = td,   
  start = list( B = rep(0, ncol(td2$PF) ),
                W = rep(0, ncol(td2$TC) ) ) )

mth.pt.u <- ulam(
  alist(
    vpm ~ dbinom(1, p),
    logit(p) <- A + PF %*% B + TC %*% W,
    A ~ dnorm(0, 3),
    B ~ dnorm(0, 5),
    W ~ dnorm(0, 5) ), 
  data = td2,   
  start = list( B = rep(0, ncol(td2$PF) ),
                W = rep(0, ncol(td2$TC) ) ), 
  chains = 4, 
  cores = 4)

N <- 1e5
post.pt <- extract.samples(mth.pt, n = N)
pred.pt <- inv_logit( 
    mean(post.pt$A) +
    td$PF %*% apply(post.pt$B, 2, mean) + 
    td$TC %*% apply(post.pt$W, 2, mean) )
pred.pt.post <- sapply(1:N, 
                       function(i) 
                         post.pt$A[i] + 
                         td$PF %*% post.pt$B[i,] + 
                         td$TC %*% post.pt$W[i,] )
pred.pt.post.p <- inv_logit(pred.pt.post)

vpm.pred.md <- apply(pred.pt.post.p, 1, median)
vpm.pred.mu <- apply(pred.pt.post.p, 1, mean)
library(devEMF)
emf(
  file = "C:\\Users\\neuro\\seadrive_root\\Emilio U\\My Libraries\\My Library\\VPM-POm spatial correlation\\VPM likelihood1.emf",
  emfPlus = TRUE,
  bg = "transparent",
  fg = "black",
  emfPlusRaster = TRUE
  )
plot(NULL, xlim = c(0,length(vpm.pred.mu)+1), ylim = c(-0.1, 1.1), ylab="VPM likelihood", 
     xlab="Cell index")
lines(1:length(vpm.pred.mu), vpm.pred.md, lwd=2, type = "l", col="black")
for (i in c(0.75,0.5,0.25,0.1,0.05)) {
  shade(apply(pred.pt.post.p, 1, PI, prob = i), 
        1:length(vpm.pred.mu), 
        col = col.alpha("black", 0.3))
}
lines(1:length(vpm.pred.mu), vpm.pred.mu, lwd=3, type="l", col=rangi2)
dev.off()

mth.p <- quap(
  alist(
    vpm ~ dbinom(1, p),
    logit(p) <- A + PF %*% B ,
    A ~ dnorm(0, 3),
    B ~ dnorm(0, 5)
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