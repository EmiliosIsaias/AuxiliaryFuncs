library(rethinking)
library(R.matlab)

fpath <- "C:\\Users\\jefe_\\seadrive_root\\Emilio U\\Für meine Gruppen\\GDrive GrohLab\\Projects\\00 Jitter\\Bayes_psth_test_1_neu.mat"

dat <- R.matlab::readMat(fpath)

d <- list(
  y = as.integer( dat$psth.u ),
  trig_flag = as.logical( dat$trigger ),
  condition = as.integer( dat$cond.id ),
  tid = as.integer(dat$trial.id)
)

mJ.1 <- ulam(
  alist(
    y ~ dzipois( s , lambda ),
    logit(s) <- alpha_s[]
  ) )