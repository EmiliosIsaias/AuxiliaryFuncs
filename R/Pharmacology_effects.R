library(rethinking)
library(R.matlab)

fpath <- "C:\\Users\\neuro\\seadrive_root\\Emilio U\\Shared with groups\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 1\\Matlab figures\\Data\\PTX_MUS_4_Stan.mat"

dat <- R.matlab::readMat(fpath)

bi_centre <- dat$bi_centre[1,]
bi_scale <- dat$bi_scale[1,]

Nt <- max( dat$tid )

d <- list(
  bi = dat$bi[,1],
  tid = as.integer( dat$tid ),
  Nt = Nt,
  mouse_id = as.integer( dat$mouse.id[,1] )
)

mPharma <- ulam(
  alist(
    treat ~ normal( mu, sigma ),
    ctr ~ normal( mu_c, sigma_c ),
    mu <- mu_c + g[treat_id],
    
    # fixed priors
    mu_c ~ normal(0, 1),
    g[treat_id] ~ normal(0, 1),
    sigma ~ dexp(1),
    sigma_c ~dexp(1)
  ) , data=d , chains=4 , cores=4 , log_lik=TRUE )

post <- link(mPharma)
params <- extract.samples( mPharma )

R.matlab::writeMat("C:\\Users\\neuro\\seadrive_root\\Emilio U\\Shared with groups\\GDrive GrohLab\\Projects\\00 SC\\SC Behaviour\\Figures\\Figure 1\\Matlab figures\\Data\\PTX and muscimol Bayes.mat", 
                   post = post, params = params, bi_centre = bi_centre, bi_scale = bi_scale)