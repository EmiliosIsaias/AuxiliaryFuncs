library(rethinking)
library(R.matlab)

fpath <- "Z:\\Jesus\\Jittering\\FULLYCurated\\05_190701_Jesus_Emili_Jittering_3800_1600_1500 yes VPM good\\KS2_newChanMap\\PSTHu_PSTHc.mat"

dat <- R.matlab::readMat(fpath)

d <- list(
  Hy = matrix( as.integer(dat$psth.c), 
               nrow = nrow(dat$psth.c), 
               ncol = ncol(dat$psth.c)),
  Nt = as.integer( nrow(dat$psth.c) ),
  U = as.integer( dim(dat$psth.u)[3] ),
  Hx =   dat$psth.u
  )

mRes <- ulam(
  alist(
    Hy ~ dzipois( f, lambda ),
    logit(f) <- t(Hx[,,u]) %*% alpha[u],
    log(lambda) <- Hx[,,u] %*% beta[u],
    
    alpha ~ dnorm( -1 , 1.5 ),
    beta ~ dnorm( 1 , 0.5 )
  ), data = d )

