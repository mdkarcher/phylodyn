int_invexp_interval <- function(t1, t2, beta0, beta1)
{
  if (beta1 == 0)
    result <- exp(-beta0) * (t2 - t1)
  else
    result <- exp(-beta0) * (exp(-beta1 * t1) - exp(-beta1 * t2)) / beta1
  return(result)
}

exp_coal_lik_init <- function(samp_times, n_sampled, coal_times)
{
  args <- gen_INLA_args(samp_times, n_sampled, coal_times)
  
  t1s <- utils::head(args$s, -1)
  t2s <- utils::tail(args$s, -1)
  
  return(list(coal_times = coal_times, t1s = t1s, t2s = t2s,
              coal_factor = args$coal_factor))
}

exp_coal_loglik <- function(betas, init)
{
  n <- length(init$coal_times)
  loglik <- -n*betas[1] + sum(betas[2] * init$coal_times)
  
  for (i in 1:length(init$t1s))
  {
    loglik <- loglik - init$coal_factor[i] * int_invexp_interval(t1 = init$t1s[i], t2 = init$t2s[i], beta0 = betas[1], beta1 = -betas[2])
  }
  
  return(loglik)
}

#' Maximum likelihood exponential model
#' 
#' @param samp_times vector of sampling times.
#' @param n_sampled vector of number sampled at \code{samp_times}.
#' @param coal_times vector of coalescent times.
#' @param betas starting values for optimization.
#'
#' @export
ml_exp_coal <- function(samp_times, n_sampled, coal_times, betas = c(0, 0))
{
  ini <- exp_coal_lik_init(samp_times, n_sampled, coal_times)
  opt <- stats::optim(betas, fn = exp_coal_loglik, init = ini, control=list(fnscale=-1))
  return(opt$par)
}