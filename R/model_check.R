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

samp_to_chisq = function(samp_times, n_sampled, logpop, grid, beta0, beta1, maxt=NULL,
                         covariates = NULL, cov_betas = NULL,
                         pow_covariates = NULL, pow_cov_betas = NULL)
  # f, beta0=1, beta1=1, maxt=Inf, pts_per_bin=1)
{
  y = exp(log_sampling_intensity(logpop = logpop, beta0 = beta0, beta1 = beta1,
                                 covariates = covariates, cov_betas = cov_betas,
                                 pow_covariates = pow_covariates,
                                 pow_cov_betas = pow_cov_betas))
  
  if (is.null(maxt))
    maxt = max(samp_times)
  maxt = min(maxt, max(grid))
  s = rep(samp_times, n_sampled)
  s = s[s <= maxt]
  grid = c(grid[grid < maxt], maxt)
  # if (maxt > max(grid))
  #  grid = c(grid, maxt)
  
  diffs = diff(grid)
  df = length(diffs)
  ys = y[1:df]
  
  expected = diffs * ys
  observed = as.numeric(table(cut(x = s, breaks = grid, right = FALSE, include.lowest = TRUE)))
  
  chisq = sum((observed - expected)^2/expected)
  
  return(list(chisq=chisq, df=df, observed=observed, expected=expected))
}


#### Single step check functions ####

samp_check_exp = function(samp_times, n_sampled, logpop, grid, betas)
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


samp_check_exp2 = function(samp_times, n_sampled, logpop, grid, betas,
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

#' Predictive check for coalescent model
#'
#' @param coal_times numeric vector of coalescent times.
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled integer vector of lineages sampled at each sampling time.
#' @param logpop numeric vector of log-effective population sizes.
#' @param grid numeric grid of time points.
#'
#' @return named numeric vector: "obs" contains the KS discrepancy for the
#'   supplied coalescent times, "rep" contains å KS discrepancy for the
#'   replicated coalescent times.
#' @export
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

#' Predictive check for sampling model
#'
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled integer vector of lineages sampled at each sampling time.
#' @param logpop numeric vector of log-effective population sizes.
#' @param grid numeric grid of time points.
#' @param betas numeric vector of log-linear coefficients for sampling model.
#' @param covariates list of vectors representing covariate values.
#' @param cov_betas numeric coefficients for each element of covariates.
#' @param pow_covariates list of vectors representing interaction covariate
#'   values.
#' @param pow_cov_betas numeric coefficients for each element of pow_covariates.
#'
#' @return named numeric vector: "obs" contains the chi-square discrepancy for
#'   the supplied sampling times, "rep" contains å chi-square discrepancy for
#'   the replicated sampling times.
#' @export
samp_check = function(samp_times, n_sampled, logpop, grid, betas,
                      covariates = NULL, cov_betas = NULL,
                      pow_covariates = NULL, pow_cov_betas = NULL)
{
  minsamp = min(samp_times)
  maxsamp = max(samp_times)
  chisq_obs = samp_to_chisq(samp_times = samp_times, n_sampled = n_sampled,
                            logpop = logpop, grid = grid,
                            beta0 = betas[1], beta1 = betas[2],
                            covariates = covariates, cov_betas = cov_betas,
                            pow_covariates = pow_covariates, pow_cov_betas = pow_cov_betas)
  
  log_int = log_sampling_intensity(logpop = logpop, beta0 = betas[1], beta1 = betas[2],
                                   covariates = covariates, cov_betas = cov_betas,
                                   pow_covariates = pow_covariates, pow_cov_betas = pow_cov_betas)
  
  int_fn = stats::stepfun(x = grid, y = exp(c(log_int[1], log_int, utils::tail(log_int, 1))))
  samp_rep = c(minsamp, pref_sample(f = int_fn, lim = c(minsamp, maxsamp)), maxsamp)
  
  chisq_rep = samp_to_chisq(samp_times = samp_rep, n_sampled = rep(1, length(samp_rep)),
                            logpop = logpop, grid = grid,
                            beta0 = betas[1], beta1 = betas[2],
                            covariates = covariates, cov_betas = cov_betas,
                            pow_covariates = pow_covariates, pow_cov_betas = pow_cov_betas)
  
  result = c(chisq_obs$chisq, chisq_rep$chisq)
  names(result) = c("obs", "rep")
  
  return(result)
}

#' Posterior predictive check for coalescent model
#'
#' @param coal_times numeric vector (for fixed-tree) or numeric vector (for
#'   inferred-tree) of coalescent times.
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled integer vector of lineages sampled at each sampling time.
#' @param logpop matrix posterior sample of log-effective population sizes.
#' @param grid numeric grid of time points.
#' @param cap integer maximum number of posterior samples to use.
#'
#' @return data.frame with two columns: "obs" contains the KS discrepancies for
#'   the supplied coalescent times, "rep" contains the KS discrepancies for the
#'   replicated coalescent times.
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
    if (is.matrix(coal_times))
      coal_vec = coal_times[i,]
    else
      coal_vec = coal_times
    
    result[i - first + 1, ] = coal_check(coal_times = coal_vec, samp_times = samp_times,
                                         n_sampled = n_sampled, logpop = as.numeric(logpop[i,]),
                                         grid = grid)
  }
  
  return(result)
}

#' Posterior predictive check for inhomogenous Poisson process sampling model
#'
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled integer vector of lineages sampled at each sampling time.
#' @param logpop matrix posterior sample of log-effective population sizes.
#' @param grid numeric grid of time points.
#' @param betas matrix posterior sample of log-linear coefficients for sampling
#'   model.
#' @param cap integer maximum number of posterior samples to use.
#' @param covariates list of vectors representing covariate values.
#' @param cov_betas numeric coefficients for each element of covariates.
#' @param pow_covariates list of vectors representing interaction covariate
#'   values.
#' @param pow_cov_betas numeric coefficients for each element of pow_covariates.
#'
#' @return data.frame with two columns: "obs" contains the chi-square
#'   discrepancies for the supplied sampling times, "rep" contains the
#'   chi-square discrepancies for the replicated sampling times.
#' @export
posterior_samp_check = function(samp_times, n_sampled, logpop, grid, betas, cap=NULL,
                                covariates = NULL, cov_betas = NULL,
                                pow_covariates = NULL, pow_cov_betas = NULL)
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
    result[i - first + 1, ] = samp_check(samp_times = samp_times,
                                         n_sampled = n_sampled, logpop = as.numeric(logpop[i,]),
                                         grid = grid, betas = as.numeric(betas[i,]),
                                         covariates = covariates, cov_betas = cov_betas,
                                         pow_covariates = pow_covariates, pow_cov_betas = pow_cov_betas)
  }
  
  return(result)
}
