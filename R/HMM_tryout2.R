library(rethinking)
library(R.matlab)

stan_model_file <-
  "C:/Users/neuro/Documents/MATLAB/AuxiliaryFuncs/Stan/HMM_2states.stan"

fpath <- "Z:\\Jesus\\Jittering\\FULLYCurated\\05_190701_Jesus_Emili_Jittering_3800_1600_1500 yes VPM good\\KS2_newChanMap\\LFP_4_R.mat"

dat <- R.matlab::readMat(fpath)
d <- list(
  N = as.integer(dat$N),
  y = drop(dat$lfp.z)
)

mLFP.1 <- stan( stan_model_file, cores=6, chains=4, data=d )
