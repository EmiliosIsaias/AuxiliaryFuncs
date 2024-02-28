library(rethinking)
library(readxl)
library(R.matlab)

fpath <- 'C:\\Users\\neuro\\seadrive_root\\Emilio U\\Shared with groups\\GDrive GrohLab\\Projects\\00 SC\\SC Anatomy\\Anatomy GADcre_VglutCre\\Analysis_raw_processed_SC_Anatomy_GAD_Vglut_Ayu.xlsx'

dat <- as.data.frame( read_xlsx( 
  path = fpath,
  sheet = 'Sheet2', 
  col_names = TRUE, 
  progress = readxl_progress() ) )

N <- 12
d <- list(
  N = N,
  O = as.integer(dat[1:N, "Overlap"]),
  Bo = as.integer(dat[1:N, "BFP"]),
  Go = as.integer(dat[1:N, "GFP"]),
  log_C = log( dat[1:N, "Total no. neurons"] ),
  gid = ifelse(dat[1:N, "Genotype"]=="Gad", 2, 1),
  C = dat[1:N, "Total no. neurons"]
)

dp <- list(
  N = N,
  C = dat$`Total no. neurons`[1:N],
  #C = dat$BFP[1:N] - dat$Overlap[1:N] + dat$GFP[1:N],
  Op = dat$Overlap[1:N]/C,
  Bop = dat$BFP[1:N]/C,
  Gop = dat$GFP[1:N]/C
)

mGlu.1 <- ulam(
  alist(
    
    Bo ~ dpois( lambda_b ),
    log( lambda_b ) <- B + b*O,
    
    Go ~ dpois( lambda_g ),
    log( lambda_g ) <- G + g*O,
    
    c(b,g) ~ dnorm( 0, 1),
    
    B ~ dnorm( mu_b, sigma_b ),
    G ~ dnorm( mu_g, sigma_g ),
    
    c(mu_b, mu_g) ~ dnorm( 0, 1),
    c(sigma_b, sigma_g) ~ dexp(1)
                               
  ), data = d, cores = 4, chains = 4 )

mGlu.2 <- ulam(
  alist(
    
    Bo ~ dpois( lambda_b ),
    log( lambda_b ) <- B + b*O,
    
    Go ~ dpois( lambda_g ),
    log( lambda_g ) <- G + g*O,
    
    C ~ dpois( lambda_c ),
    log( lambda_c ) <- B + G + (b+g)*O,
    
    B ~ dnorm( 15, 5),
    G ~ dnorm( 60, 20),
    c(b,g) ~ dnorm( 0, 1)
    
    
  ), data = d, cores = 4, chains = 4 )

mGlu.3 <- ulam(
  alist(
    
    Bo ~ dpois( lambda_b ),
    log( lambda_b ) <- log_C - G[gid] - ob[gid]*O,
    
    Go ~ dpois( lambda_g ),
    log( lambda_g ) <- log_C - B[gid] - og[gid]*O,
    
    O ~ dpois( lambda_o ),
    log( lambda_o ) <- log_C - B[gid] - G[gid],
    
    B[gid] ~ dnorm( 0, 2 ),
    G[gid] ~ dnorm( 0, 3 ),
    ob[gid] ~ dnorm( 0, 1 ),
    og[gid] ~ dnorm( 0, 1 )
    
  ), data = d, chains = 4, cores = 4 )

precis( mGlu.3 , depth = 2 )

post.3 <- link( mGlu.3 )
params.3 <- extract.samples( mGlu.3 )

matrix( c( d$Bo, apply( post.3$lambda_b, 2, mean) ), ncol = 2)
matrix( c( d$Go, apply( post.3$lambda_g, 2, mean) ), ncol = 2)

plot( d$Bo, apply( post.3$lambda_b, 2, mean))
abline(a = 0, b = 1, col = "black", lwd = 2)
for (i in 1:N) {
  lines( c(d$Bo[i], d$Bo[i]), PI( post.3$lambda_b[,i], prob = 0.9) )
}

plot( d$Go, apply( post.3$lambda_g, 2, mean))
abline(a = 0, b = 1, col = "black", lwd = 2)
for (i in 1:N) {
  lines( c(d$Go[i], d$Go[i]), PI( post.3$lambda_g[,i], prob = 0.9) )
}

plot( d$O, apply( post.3$lambda_o, 2, mean) )
abline(a = 0, b = 1, col = "black", lwd = 2)
for (i in 1:N) {
  lines( c(d$O[i], d$O[i]), PI( post.3$lambda_o[,i], prob = 0.9) )
}

matrix(c(d$Bo, 
         exp(d$log_C - mean( params.3$G )), 
         exp(d$log_C - mean( params.3$G ) - mean( params.3$b)*d$O),
         exp(d$log_C) ), ncol = 4)

matrix(c(d$Go, 
         exp(d$log_C - mean( params.3$B )), 
         exp(d$log_C - mean( params.3$B ) - mean( params.3$g)*d$O),
         exp(d$log_C) ), ncol = 4)

# Overlap estimation given the other numbers
mGlu.4 <- quap(
  alist(    
    Bo ~ dpois( lambda_b ),
    lambda_b <- C - g*Go,
    
    Go ~ dpois( lambda_g ),
    lambda_g <- C - b*Bo,
    
    C ~ dpois( lambda_c ),
    lambda_c <- b*Bo + g*Go + O,
    
    O ~ dpois( lambda_o ),
    lambda_o <- b*Bo,
    
    c(b, g) ~ dlnorm( 0, 1 )
    
  ), data = d, start = list(b=1, g=1))

precis( mGlu.4 )

# Winner of this contest.
mGlu.6 <- ulam(
  alist(
    
    O ~ dpois( lambda_o ),
    lambda_o <- f[gid]*Bo,
    
    C ~ dpois( lambda_c ),
    lambda_c <- B[gid] + Go,
    
    B[gid] ~ dlnorm( 0, 2 ),
    f[gid] ~ dlnorm( 0.5, 2 )
    
  ), data = d, chains = 4, cores = 4 )

precis( mGlu.6, depth = 2 )

post.6 <- link( mGlu.6 )
params.6 <- extract.samples( mGlu.6 )

