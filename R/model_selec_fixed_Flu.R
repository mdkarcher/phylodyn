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
  
  ## model 1: no pref sampling (conditional)
  fit1 <- BNPR_fixed(data = phy, lengthout = 100)
  ## model 2: pref sampling (no covariates)
  fit2 <- BNPR_PS(data = phy, lengthout = 100)
  ## model 3: pref sampling, time covariate
  fit3 <- BNPR_PS(data = phy,  fns = list(tcov), log_fns = FALSE, lengthout = 100)
  ## model 4: pref sampling with SEASONAL covariates
  fit4 <- BNPR_PS(data = phy, lengthout = 100, 
                  fns = list(winter_ind, autumn_ind, summer_ind), log_fns = FALSE)
  ## model 5: conditional model + nonparametric model of the sampling rate
  fit5 <- BNPR(data = phy, lengthout = 100)
  aux.fit5 <- BNPR_samp_only(data = phy, lengthout = 100)
  ##
  mal1 <- get_marginal_likelihood(fit1)
  mal2 <- get_marginal_likelihood(fit2)
  mal3 <- get_marginal_likelihood(fit3)
  mal4 <- get_marginal_likelihood(fit4)
  logadjustment <- as.numeric(aux.fit5$internals$result$mlik)
  mal5 <- get_marginal_likelihood(fit5, adj = logadjustment)
  MaLs <- data.frame(method = mal1$method,
                     m1 = mal1$marginal_likelihood,
                     m2 = mal2$marginal_likelihood,
                     m3 = mal3$marginal_likelihood,
                     m4 = mal4$marginal_likelihood,
                     m5 = mal5$marginal_likelihood)
  return(
    list(
      log_marginal_likelihoods = MaLs,
      beta2_estimates = get_beta2_estimates(fit3)
    )
  )
}

############ Data and fitting
lengthout <- 100
flu <- read.nexus("data/USACanada_DNA_skygrid_max_clade.nexus")
gene.flu <- summarize_phylo(flu)

width <- max(gene.flu$coal_times)/lengthout
epochs <- rep(width, lengthout + 1)
midpts <- cumsum(epochs) - epochs/2

## covariate functions
winter_ind  <- function(t) {
  as.numeric(t %% 1.0 < 0.25)
}

autumn_ind <- function(t) {
  as.numeric((t %% 1.0 >= 0.25) & (t %% 1.0 < 0.50))
}

summer_ind <- function(t) {
  as.numeric((t %% 1.0 >= 0.50) & (t %% 1.0 < 0.75))
}

tcov <- function(t) -t

### Results

fit_and_report(flu) ## best model is seasonal model

