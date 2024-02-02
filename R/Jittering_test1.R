library(rethinking)
library(R.matlab)

fpath <- "C:\\Users\\jefe_\\seadrive_root\\Emilio U\\Für meine Gruppen\\GDrive GrohLab\\Projects\\00 Jitter\\Bayes_psth_test.mat"
dat <- R.matlab::readMat(fpath)

d <- list(
  y = matrix(as.integer(dat$psth.stack), nrow = 6560, ncol = 4500),
  cid = as.integer(dat$cl.id),
  tid = as.integer(dat$tr.id),
  tm = 1:4500
)

m1 <- ulam(
  alist(
    y ~ dzipois( lambda_s ),
    log(lambda_s) <- al,
    al ~ normal( 1 , 2 )
  ), data = d, chains = 4, cores = 4)

# Estimating triggers and firing rates from a 0.1 ms bin PSTH per cluster per trial
m1 <- ulam(
  alist(
    y ~ dzipois( t , lambda_s ),
    logit(t) <- at[tm],
    log(lambda_s) <- al[cid,tid],
    vector[4500]:at ~ dnorm( 0 , 10 ),
    matrix[164,tid]:al ~ dnorm( 1 , 0.5 )
  ) , data = d , chains = 4, cores = 4)
