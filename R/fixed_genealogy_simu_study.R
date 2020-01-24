## Set up 
## Given a tree size (n) and value for beta2 (b2), simulate **both** (a) sampling times 's' and (b) tree conditional on s;
## Fit a model that accounts for preferential sampling (BNPR1) and a model that accounts for PS **and** includes the right covariate (BNPR2)
## Extract the posterior mean and CI for beta2 from BNPR2 and compute the log Bayes factor between BNPR1 and BNPR2

Ncores <- 20

library(phylodyn)
source("BNPR_PS_fixedBeta.R")
# INLA:::inla.dynload.workaround()

### Prelim stuff
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
#
generate_data <- function(n, scenario, beta1, beta2, samp_end = 60){
  
  traj1 <- function(t) logistic_traj(t, offset = 6) / 5
  traj2 <- function(t) exp_traj(t = t, scale = 1, rate = 0.05)
  
  ## scenario 1: no pref sampling
  s1 <- function(){
    trajc <- function(t) 1 * t^0
    trajc <- Vectorize(trajc)
    Cprop.unif  <- n/integrate(traj_beta, 0, samp_end, traj = trajc, beta = 1)$value
    uniform.times <- pref_sample(trajc, c = Cprop.unif, lim = c(0, samp_end), beta = 1)
    return(coalsim(samp_times = uniform.times,
                   n_sampled = rep(1, length(uniform.times)), traj = traj1)) 
  }
  ## scenario 2: pref sampling, no covariates
  s2 <-  function(){
    trajc <- function(t) traj1(t)^beta1
    trajc <- Vectorize(trajc)
    Cprop.pref  <- n/integrate(traj_beta, 0, samp_end, traj = traj1, beta = 1)$value
    pref.times <- pref_sample(traj1, c = Cprop.pref, lim = c(0, samp_end), beta = 1)
    return(coalsim(samp_times = pref.times,
                   n_sampled = rep(1, length(pref.times)), traj = traj1))
  }
  ## scenario 3: pref sampling with covariates
  s3 <- function(){
    trajc <- function(t) traj1(t)^beta1 * traj2(t)^beta2
    trajc <- Vectorize(trajc)
    Cprop.prefCov  = n/integrate(traj_beta, 0, samp_end, traj = trajc, beta = 1)$value
    prefCov.times = pref_sample(trajc, c = Cprop.prefCov, lim = c(0, samp_end), beta = 1)
    return(coalsim(samp_times = prefCov.times,
                   n_sampled = rep(1, length(prefCov.times) ), traj = traj1))
  }
  ## scenario 4: unrelated [to both Net(t) and the covariate(s)] pref sampling
  s4 <- function(){
    trajc <- function(t) 10* abs(cos(.05 * t)) ## trough in the middle of the sampling span
    trajc <- Vectorize(trajc)
    Cprop.Unrelated  = n/integrate(traj_beta, 0, samp_end, traj = trajc, beta = 1)$value
    prefUnrelated.times = pref_sample(trajc, c = Cprop.Unrelated, lim = c(0, samp_end), beta = 1)
    return(coalsim(samp_times = prefUnrelated.times,
                   n_sampled = rep(1, length(prefUnrelated.times)), traj = traj1))
  }
  
  out <- switch(scenario,
                "1" = s1(),
                "2" = s2(),
                "3" = s3(),
                "4" = s4())
  return(out)
}
#
fit_and_report <- function(phy){
  
  traj1 <- function(t) logistic_traj(t, offset = 6) / 5
  traj2 <- function(t) exp_traj(t = t, scale = 1, rate = 0.05)
  
  ## model 1: no pref sampling (conditional)
  fit1 <- BNPR_fixed(data = phy, lengthout = 100)
  ## model 2: pref sampling, no covariates
  fit2 <- BNPR_PS(data = phy, lengthout = 100)
  ## model 3: pref sampling with covariates
  fit3 <- BNPR_PS(data = phy, lengthout = 100, fns = list(traj2))
  ## model 4: conditional model + nonparametric model of the sampling rate
  fit4 <- BNPR(data = phy, lengthout = 100)
  aux.fit4 <- BNPR_samp_only(data = phy, lengthout = 100)
  ##
  mal1 <- get_marginal_likelihood(fit1)
  mal2 <- get_marginal_likelihood(fit2)
  mal3 <- get_marginal_likelihood(fit3)
  logadjustment <- as.numeric(aux.fit4$internals$result$mlik)
  mal4 <- get_marginal_likelihood(fit4, adj = logadjustment)
  MaLs <- data.frame(method = mal1$method,
                     m1 = mal1$marginal_likelihood,
                     m2 = mal2$marginal_likelihood,
                     m3 = mal3$marginal_likelihood,
                     m4 = mal4$marginal_likelihood)
  measures <- lapply(
    list(model1 = fit1, model2 = fit2, model3 = fit3, model4 = fit4),
    function(fit){
      return(data.frame(
        MAD = mean(abs(fit$effpop - traj1(fit$x))),
        MCIW = mean(fit$effpop975- fit$effpop025)
      ))
    }
  )
  measures.df <- do.call(rbind, measures)
  measures.df$model <- rownames(measures.df)
  rownames(measures.df) <- NULL
  return(
    list(
      error_measures = measures.df,
      log_marginal_likelihoods = MaLs,
      beta2_estimates = get_beta2_estimates(fit3)
    )
  )
}
###
## Now simulation-specific function
simulate_fit_and_report <- function(n, beta2, beta1 = 1, samp_end = 60){
  
  s1.data <- generate_data(n = n, scenario = 1,
                           beta1 = beta1, beta2 = beta2, samp_end = samp_end)
  s2.data <- generate_data(n = n, scenario = 2,
                           beta1 = beta1, beta2 = beta2, samp_end = samp_end)
  s3.data <- generate_data(n = n, scenario = 3,
                           beta1 = beta1, beta2 = beta2, samp_end = samp_end)
  s4.data <- generate_data(n = n, scenario = 4,
                           beta1 = beta1, beta2 = beta2, samp_end = samp_end)
  
  ## Now fit the models and extract relevant info
  stamp_scenario <- function(x, scenario_number){
    out <- lapply(x, function(sub.x) data.frame(sub.x, scenario = paste0(scenario_number)))
    return(out)
  }

  results.scenario.1 <- stamp_scenario (fit_and_report(s1.data), 1)
  results.scenario.2 <- stamp_scenario (fit_and_report(s2.data), 2)
  results.scenario.3 <- stamp_scenario (fit_and_report(s3.data), 3)
  results.scenario.4 <- stamp_scenario (fit_and_report(s4.data), 4)
  
  error.measures <- do.call(rbind,
                            list(results.scenario.1$error_measures,
                                 results.scenario.2$error_measures,
                                 results.scenario.3$error_measures,
                                 results.scenario.4$error_measures))
  error.measures <- data.frame(sample_size = n, true_beta2 = beta2, error.measures)
  #
  log.marginal.likelihoods <- do.call(rbind,
                                      list(results.scenario.1$log_marginal_likelihoods,
                                           results.scenario.2$log_marginal_likelihoods,
                                           results.scenario.3$log_marginal_likelihoods,
                                           results.scenario.4$log_marginal_likelihoods))
  log.marginal.likelihoods <- data.frame(sample_size = n, true_beta2 = beta2,
                                         log.marginal.likelihoods)
  #
  beta2.estimates <- do.call(rbind,
                             list(results.scenario.1$beta2_estimates,
                                  results.scenario.2$beta2_estimates,
                                  results.scenario.3$beta2_estimates,
                                  results.scenario.4$beta2_estimates
                             ))
  beta2.estimates <- data.frame(sample_size = n, true_beta2 = beta2, beta2.estimates)
  return(
    list(
      error_measures =  error.measures,
      marginal_likelihoods = log.marginal.likelihoods,
      beta2_estimates =  beta2.estimates
    )
  )
}
#
rep_simulationfit <- function(Ns, b2, nrep, cores = 20){
  return(
    do.call(
      rbind,
      parallel::mclapply(1:nrep, function(i){
        res <- simulate_fit_and_report(n = Ns, beta2 = b2)
        return(
          lapply(res, function(x) data.frame(x, replicate = i))
        )
      }, mc.cores = cores)
    )
  )
}

######
simu.grid <- expand.grid(
  sample_size = c(50, 100, 500, 1000),
  tbeta2 = c(.5, 1, 5, 10)
)

K <- nrow(simu.grid)
Error <- vector(K, mode = "list")
lMaL <- vector(K, mode = "list")
Beta2 <- vector(K, mode = "list")
Nrep <- 500

system("echo '' > simu_monitor.txt")

simu.time <- system.time(
  for(k in 1:K){
    write(paste("Doing simulation ", k, " of ", K, "\n", sep = ""),
          file = "simu_monitor.txt", append = TRUE)
    Simu <- rep_simulationfit(Ns = simu.grid[k, ]$sample_size,
                              b2 = simu.grid[k, ]$tbeta2,
                              nrep = Nrep, cores = Ncores)
    Error[[k]] <-  do.call(rbind, Simu[, 1])
    lMaL[[k]] <- do.call(rbind, Simu[, 2])
    Beta2[[k]] <- do.call(rbind, Simu[, 3])
  }
)
simu.time

filter_problematic <- function(dt){
  dt[-grep("integrate", dt$sample_size), ]
}

ErrorMeasures.dt <- do.call(rbind, Error)
# ErrorMeasures.dt <- filter_problematic(ErrorMeasures.dt)

ErrorMeasures.dt

lMaL.dt <- do.call(rbind, lMaL)
# lMaL.dt <- filter_problematic(lMaL.dt)
lMaL.dt

Beta2.dt <- do.call(rbind, Beta2)
# Beta2.dt <- filter_problematic(Beta2.dt)
Beta2.dt

write.csv(ErrorMeasures.dt, "simu_study_fixed_genealogy_MAD_MCIW.csv", row.names = FALSE)
write.csv(lMaL.dt, "simu_study_fixed_genealogy_logMarginalLikelihoods.csv", row.names = FALSE)
write.csv(Beta2.dt, "simu_study_fixed_genealogy_Beta2_estimates.csv", row.names = FALSE)