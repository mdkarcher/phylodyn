library(ape)
library(phylodyn)
source("BNPR_PS_fixedBeta.R")

############ Functions
get_beta2_estimates <- function(fit){
  dt <-  fit$result$summary.fixed[2, ][c(1, 3, 5)]
  colnames(dt) <- c("mean", "lwr", "upr")
  rownames(dt) <- NULL
  return(dt)
}
#
get_marginal_likelihood <- function(fit, adj = 0){
  raw <- summary(fit$result)$mlik
  return(data.frame(
    method = c("Integration", "Gaussian"),
    marginal_likelihood = as.numeric(raw) + adj
  ))
}
##
fit_and_report <- function(phy){
  traj2 <- function(t) -t
  traj3 <- function(t) -t*t
  ## model 1: no pref sampling (conditional)
  fit1 <- BNPR_fixed(data = phy, lengthout = 100)
  ## model 2: pref sampling, no covariates
  fit2 <- BNPR_PS(data = phy, lengthout = 100)
  ## model 3a: pref sampling with -t
  fit3a <- BNPR_PS(data = phy, lengthout = 100, fns = list(traj2), log_fns = FALSE)
  ## model 3c: pref sampling with -t^2
  fit3b <- BNPR_PS(data = phy, lengthout = 100, fns = list(traj3), log_fns = FALSE)
  ## model 3c: pref sampling with both covariates
  fit3c <- BNPR_PS(data = phy, lengthout = 100, fns = list(traj2, traj3), log_fns = FALSE)
  ## model 4: conditional model + nonparametric model of the sampling rate
  fit4 <- BNPR(data = phy, lengthout = 100)
  aux.fit4 <- BNPR_samp_only(data = phy, lengthout = 100)
  ##
  mal1 <- get_marginal_likelihood(fit1)
  mal2 <- get_marginal_likelihood(fit2)
  mal3a <- get_marginal_likelihood(fit3a)
  mal3b <- get_marginal_likelihood(fit3b)
  mal3c <- get_marginal_likelihood(fit3c)
  logadjustment <- as.numeric(aux.fit4$internals$result$mlik)
  mal4 <- get_marginal_likelihood(fit4, adj = logadjustment)
  MaLs <- data.frame(method = mal1$method,
                     m1 = mal1$marginal_likelihood,
                     m2 = mal2$marginal_likelihood,
                     m3a = mal3a$marginal_likelihood,
                     m3b = mal3b$marginal_likelihood,
                     m3c = mal3c$marginal_likelihood,
                     m4 = mal4$marginal_likelihood)
  return(
    list(
      log_marginal_likelihoods = MaLs
      # beta2_estimates = get_beta2_estimates(fit3)
    )
  )
}

############ Data and fitting
##  Trees have polytomies, which we resolve randomly.
## Analyses were run many times to ensure no effect of this step is negligible

poly.SLE <- read.tree("data/prefsample_200tt_dated.tree") ## Sierra Leone
resolved.SLE <- multi2di(poly.SLE, random = TRUE)
fit_and_report(resolved.SLE)

poly.LBR <- read.tree("data/prefsample_Liberiatt_dated.tree") ## Liberia
resolved.LBR <- multi2di(poly.LBR, random = TRUE)
fit_and_report(resolved.LBR)

