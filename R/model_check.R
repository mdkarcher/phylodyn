#### Helper functions ####

integrate_step_fun = function(fun, from, to)
{
  #fun = stats::stepfun(x = breaks, y = c(min(y), y, max(y)))
  breaks = stats::knots(fun)
  #levels = environment(fun)$y
  
  mask = breaks > from & breaks < to
  pts = c(from, breaks[mask], to)
  diffs = diff(pts)
  midpts = pts[-1] - diffs/2
  
  segments = diffs * fun(midpts)
  
  return(sum(segments))
}

#### Conversion functions ####

coal_to_exp = function(coal_times, samp_times, n_sampled, logpop, grid, maxt=Inf)
{
  coal_times = coal_times[coal_times <= maxt]
  y = exp(-logpop)
  Fn = stats::stepfun(x = grid, y = c(min(y), y, max(y)))
  
  # first case special: first sample to first coalescent
  sampling_mask = samp_times < coal_times[1] & samp_times > samp_times[1]
  knots = c(samp_times[1], samp_times[sampling_mask], coal_times[1])
  e = 0
  for (j in 2:length(knots))
  {
    k = sum(n_sampled[samp_times < knots[j]])
    c = k * (k-1) / 2
    e = e + c * integrate_step_fun(fun = Fn, from = knots[j-1], to = knots[j])
    #e = e + c * integrate(invf, lower=knots[j-1], upper=knots[j])$value
  }
  result = e
  
  # 2nd onwards normal: (i-1)th coalescent to (i)th coalescent
  i=2
  while (i <= length(coal_times))
  {
    sampling_mask = samp_times < coal_times[i] & samp_times > coal_times[i-1]
    knots = c(coal_times[i-1], samp_times[sampling_mask], coal_times[i])
    e = 0
    for (j in 2:length(knots))
    {
      k = sum(n_sampled[samp_times < knots[j]]) - sum(coal_times < knots[j])
      c = k * (k-1) / 2
      e = e + c * integrate_step_fun(fun = Fn, from = knots[j-1], to = knots[j])
    }
    result = c(result, e)
    i = i + 1
  }
  return(result)
}

samp_to_exp = function(samp_times, n_sampled, logpop, grid, beta0=0, beta1=1, maxt=Inf)
{
  y = exp(beta0 + beta1 * logpop)
  Fn = stats::stepfun(x = grid, y = c(min(y), y, max(y)))
  
  #f_beta = function(x) return(exp(beta0) * f(x)^beta1)
  s = rep(samp_times, n_sampled)
  s = s[s <= maxt]
  
  # first case special: time 0 to first sample
  if (s[1] > 0)
    result = integrate_step_fun(fun = Fn, from=0, to=s[1])
  else
    result = NULL
  
  # 2nd onwards normal: (i-1)th sample to (i)th sample
  i=2
  while (i <= length(s))
  {
    result = c(result, integrate_step_fun(fun = Fn, from=s[i-1], to=s[i]))
    i = i + 1
  }
  return(result)
}

exp_to_KS = function(exps, rm_zero=FALSE)
{
  if (rm_zero)
    exps = exps[exps != 0]
  return(stats::ks.test(exps, "pexp"))
}

#### Single step check functions ####

coal_check = function(coal_times, samp_times, n_sampled, logpop, grid)
{
  exps_obs = coal_to_exp(coal_times = coal_times, samp_times = samp_times,
                         n_sampled = n_sampled, logpop = logpop,
                         grid = grid)
  KSobs = exp_to_KS(exps_obs)$statistic
  
  traj = stats::stepfun(x = grid, y = exp(c(logpop[1], logpop, utils::tail(logpop,1))))
  coal_rep = coalsim(samp_times = samp_times, n_sampled = n_sampled, traj = traj)
  exps_rep = coal_to_exp(coal_times = coal_rep$coal_times, samp_times = samp_times,
                         n_sampled = n_sampled, logpop = logpop,
                         grid = grid)
  KSrep = exp_to_KS(exps_rep)$statistic
  
  result = c(KSobs, KSrep)
  names(result) = c("obs", "rep")
  
  return(result)
}

samp_check = function(samp_times, n_sampled, logpop, grid, betas)
{
  minsamp = min(samp_times)
  maxsamp = max(samp_times)
  exps_obs = samp_to_exp(samp_times = samp_times, n_sampled = n_sampled,
                         logpop = logpop, grid = grid,
                         beta0 = betas[1], beta1 = betas[2])
  KSobs = exp_to_KS(exps_obs)$statistic
  
  traj = stats::stepfun(x = grid, y = exp(c(logpop[1], logpop, utils::tail(logpop,1))))
  samp_rep = c(minsamp, pref_sample(f = traj, lim = c(minsamp, maxsamp),
                                    c = exp(betas[1]), beta = betas[2]), maxsamp)
  exps_rep = samp_to_exp(samp_times = samp_rep, n_sampled = rep(1, length(samp_rep)),
                         logpop = logpop, grid = grid,
                         beta0 = betas[1], beta1 = betas[2])
  KSrep = exp_to_KS(exps_rep)$statistic
  
  result = c(KSobs, KSrep)
  names(result) = c("obs", "rep")
  
  return(result)
}

#### Covariates version ####

log_sampling_intensity = function(logpop, beta0, beta1,
                                  covariates = NULL, cov_betas = NULL,
                                  pow_covariates = NULL, pow_cov_betas = NULL)
{
  result = beta0 + logpop * beta1
  
  if (!is.null(covariates))
  {
    for (i in 1:length(covariates))
    {
      result = result + covariates[[i]] * cov_betas[i]
    }
  }
  
  if (!is.null(pow_covariates))
  {
    for (i in 1:length(pow_covariates))
    {
      result = result + pow_covariates[[i]] * logpop * pow_cov_betas[i]
    }
  }
  
  return(result)
}

samp_to_exp2 = function(samp_times, n_sampled, logpop, grid, beta0, beta1,
                        covariates = NULL, cov_betas = NULL,
                        pow_covariates = NULL, pow_cov_betas = NULL, maxt=Inf)
{
  y = exp(log_sampling_intensity(logpop = logpop, beta0 = beta0, beta1 = beta1,
                                 covariates = covariates, cov_betas = cov_betas,
                                 pow_covariates = pow_covariates,
                                 pow_cov_betas = pow_cov_betas))
  Fn = stats::stepfun(x = grid, y = c(min(y), y, max(y)))
  
  #f_beta = function(x) return(exp(beta0) * f(x)^beta1)
  s = rep(samp_times, n_sampled)
  s = s[s <= maxt]
  
  # first case special: time 0 to first sample
  if (s[1] > 0)
    result = integrate_step_fun(fun = Fn, from=0, to=s[1])
  else
    result = NULL
  
  # 2nd onwards normal: (i-1)th sample to (i)th sample
  i=2
  while (i <= length(s))
  {
    result = c(result, integrate_step_fun(fun = Fn, from=s[i-1], to=s[i]))
    i = i + 1
  }
  return(result)
}


samp_check2 = function(samp_times, n_sampled, logpop, grid, betas,
                       covariates = NULL, cov_betas = NULL,
                       pow_covariates = NULL, pow_cov_betas = NULL)
{
  minsamp = min(samp_times)
  maxsamp = max(samp_times)
  exps_obs = samp_to_exp2(samp_times = samp_times, n_sampled = n_sampled,
                          logpop = logpop, grid = grid,
                          beta0 = betas[1], beta1 = betas[2],
                          covariates = covariates, cov_betas = cov_betas,
                          pow_covariates = pow_covariates, pow_cov_betas = pow_cov_betas)
  KSobs = exp_to_KS(exps_obs)$statistic
  
  log_int = log_sampling_intensity(logpop = logpop, beta0 = betas[1], beta1 = betas[2],
                                   covariates = covariates, cov_betas = cov_betas,
                                   pow_covariates = pow_covariates, pow_cov_betas = pow_cov_betas)
  
  int_fn = stats::stepfun(x = grid, y = exp(c(log_int[1], log_int, utils::tail(log_int, 1))))
  samp_rep = c(minsamp, pref_sample(f = int_fn, lim = c(minsamp, maxsamp)), maxsamp)
  
  exps_rep = samp_to_exp2(samp_times = samp_rep, n_sampled = rep(1, length(samp_rep)),
                          logpop = logpop, grid = grid,
                          beta0 = betas[1], beta1 = betas[2],
                          covariates = covariates, cov_betas = cov_betas,
                          pow_covariates = pow_covariates, pow_cov_betas = pow_cov_betas)
  KSrep = exp_to_KS(exps_rep)$statistic
  
  result = c(KSobs, KSrep)
  names(result) = c("obs", "rep")
  
  return(result)
}

#### User Facing Functions ####

#' Posterior predictive check for coalescent model
#' 
#' @param coal_times numeric vector of coalescent times.
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled integer vector of lineages sampled at each sampling time.
#' @param logpop matrix posterior sample of log-effective population sizes.
#' @param grid numeric grid
#' @param cap integer maximum number of posterior samples to use.
#'   
#' @export
posterior_coal_check = function(coal_times, samp_times, n_sampled, logpop, grid, cap=NULL)
{
  niter = nrow(logpop)
  result = data.frame(obs = rep(0, min(cap, niter)), rep = rep(0, min(cap, niter)))
  
  if (is.null(cap) || cap < 1)
  {
    first = 1
  }
  else
  {
    first = niter - cap + 1
  }
  
  for (i in first:niter)
  {
    result[i - first + 1, ] = coal_check(coal_times = coal_times, samp_times = samp_times,
                                         n_sampled = n_sampled, logpop = logpop[i,],
                                         grid = grid)
    # exps_obs = coal_to_exp(coal_times = coal_times, samp_times = samp_times,
    #                        n_sampled = n_sampled, logpop = logpop[i,],
    #                        grid = grid)
    # KSobs[i] = exp_to_KS(exps_obs)$statistic
    # 
    # traj = stats::stepfun(x = grid, y = exp(c(logpop[i,1], logpop[i,], utils::tail(logpop[i,],1))))
    # coal_rep = coalsim(samp_times = samp_times, n_sampled = n_sampled, traj = traj)$coal_times
    # exps_rep = coal_to_exp(coal_times = coal_rep$coal_times, samp_times = samp_times,
    #                        n_sampled = n_sampled, logpop = logpop[i,],
    #                        grid = grid)
    # KSrep[i] = exp_to_KS(exps_rep)$statistic
  }
  
  return(result)
}

#' Posterior predictive check for inhomogenous Poisson process sampling model
#' 
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled integer vector of lineages sampled at each sampling time.
#' @param logpop matrix posterior sample of log-effective population sizes.
#' @param grid numeric grid
#' @param betas matrix posterior sample of log-linear coefficients for sampling
#'   model.
#' @param cap integer maximum number of posterior samples to use.
#'   
#' @export
posterior_samp_check = function(samp_times, n_sampled, logpop, grid, betas, cap=NULL)
{
  niter = nrow(logpop)
  result = data.frame(obs = rep(0, min(cap, niter)), rep = rep(0, min(cap, niter)))
  # minsamp = min(samp_times)
  # maxsamp = max(samp_times)
  # KSobs = rep(0, min(cap, niter))
  # KSrep = rep(0, min(cap, niter))
  
  if (is.null(cap) || cap < 1)
  {
    first = 1
  }
  else
  {
    first = niter - cap + 1
  }
  
  for (i in first:niter)
  {
    result[i - first + 1, ] = samp_check(samp_times = samp_times,
                                         n_sampled = n_sampled, logpop = logpop[i,],
                                         grid = grid, betas = betas[i,])
    # exps_obs = samp_to_exp(samp_times = samp_times, n_sampled = n_sampled,
    #                        logpop = logpop[i,], grid = grid,
    #                        beta0 = betas[1], beta1 = betas[2])
    # KSobs[i] = exp_to_KS(exps_obs)$statistic
    # 
    # traj = stats::stepfun(x = grid, y = exp(c(logpop[i,1], logpop[i,], utils::tail(logpop[i,],1))))
    # samp_rep = c(minsamp, pref_sample(f = traj, lim = c(minsamp, maxsamp),
    #                                   c = exp(betas[1]), beta = betas[2]), maxsamp)
    # exps_rep = samp_to_exp(samp_times = samp_rep, n_sampled = rep(1, length(samp_rep)),
    #                        logpop = logpop[i,], grid = grid,
    #                        beta0 = betas[1], beta1 = betas[2])
    # KSrep[i] = exp_to_KS(exps_rep)$statistic
  }
  
  return(result)
}
