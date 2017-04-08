#### Elliptical slice sampler by Murray et~al (2010) ####

# inputs:
#   q_cur: initial state of the parameter
#   l_cur: initial log-likelihood
#   loglik: log-likelihood function of q
#   cholC: Cholesky decomposition (upper triangular matrix) of covariance matrix of Gaussian prior
# outputs:
#   q: new state of the parameter following N(q;0,Cov)*lik
#   u: log-likelihood of new state
#   Ind: proposal acceptance indicator

ESS = function(q_cur, l_cur, loglik, cholC, ...)
{  
  # choose ellipse
  nu = crossprod(cholC, stats::rnorm(length(q_cur)))
  
  # log-likelihood threshold
  u = stats::runif(1)
  logy <- l_cur + log(u)
  
  # draw a initial proposal, also defining a bracket
  t = 2*pi*stats::runif(1)
  t_min <- t-2*pi
  t_max <- t
  
  q <- q_cur*cos(t) + nu*sin(t)
  l <- loglik(q, ...)
  
  while (l < logy)
  {
    # shrink the bracket and try a new point
    if (t < 0)
    {
      t_min <- t
    }
    else
    {
      t_max <- t
    }
    
    t <- stats::runif(1, t_min, t_max)
    q <- q_cur*cos(t) + nu*sin(t)
    
    l <- loglik(q, ...)
  }
  
  return(list(q=q, u=l, Ind=1))
}

ESS_wrapper = function(lik_init, loglik, l_cur, f, kappa, cholC, invC, alpha, beta)
{
  #ll = function(f) coal_loglik(init = lik_init, f = f)
  res = ESS(q_cur = f, l_cur = l_cur, loglik = loglik, cholC = cholC/sqrt(kappa))
  pos_summ = list(loglik = res$u)
  pos_summ$logpri = log_mvnorm_prior(x = res$q, prec = invC * kappa) +
    log_kappa_prior(kappa = kappa, alpha = alpha, beta = beta)
  pos_summ$logpos = pos_summ$logpri + pos_summ$loglik
  res$pos_summ = pos_summ
}

ESS_old = function(q_cur, l_cur, loglik, cholC)
{  
  # choose ellipse
  nu = t(cholC) %*% stats::rnorm(length(q_cur))
  
  # log-likelihood threshold
  u = stats::runif(1)
  logy <- l_cur + log(u)
  
  # draw a initial proposal, also defining a bracket
  t = 2*pi*stats::runif(1)
  t_min <- t-2*pi
  t_max <- t
  
  while (1)
  {
    q <- q_cur*cos(t) + nu*sin(t)
    l = loglik(q)
    if (l > logy)
    {
      return(list(q=q, l=l, Ind=1))
    }
    # shrink the bracket and try a new point
    if (t < 0)
    {
      t_min <- t
    }
    else
    {
      t_max <- t
    }
    
    t = stats::runif(1, t_min, t_max)
  }
}

#### Metropolis-Adjusted Langevin (MALA) Algorithm ####
# This function generates one sample given previous state.

# inputs:
#   q_cur: initial state of the parameter
#   u_cur, du_cur: initial potential energy and its gradient
#   U:=-log(density(q)), potential function of q, or its gradient
#   eps: step size
# outputs:
#   q: new state of the parameter
#   u, du: new potential energy and its gradient
#   Ind: proposal acceptance indicator

MALA = function (q_cur, u_cur, du_cur, U, eps=.2)
{
  # initialization
  q = q_cur
  D = length(q)
  u = u_cur
  du = du_cur
  
  # sample momentum
  p = stats::rnorm(D)
  
  # calculate current energy
  E_cur = u + sum(p^2)/2
  
  # Make a half step for momentum
  p = p - eps/2 * du
  
  # Make a full step for the position
  q = q + eps * p
  
  du = U(q, grad = TRUE)$dlogpos
  # Make a half step for momentum at the end
  p = p - eps/2 * du
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  pos_summ = U(q)
  u = pos_summ$logpos
  E_prp = u + sum(p^2)/2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  logAP = -E_prp + E_cur
  
  if( is.finite(logAP) && (log(stats::runif(1)) < min(0, logAP)) ) 
    return (list(q = q, u = u, du = du, Ind = 1, pos_summ = pos_summ))
  else 
    return (list(q = q_cur, u = u_cur, du = du_cur, Ind = 0, pos_summ = U(q_cur)))
}

#### adaptive Metropolis-Adjusted Langevin (aMALA) Algorithm ####
# This is adaptive block updating GMRF by Knorr-Held and Rue (2002), equivalent to Riemannian MALA by Girolami and Calderhead (2011).
# This function generates one sample given previous state.

# inputs:
#   q_cur: initial state of the parameter
#   u_cur: initial potential energy
#   U:=-log(density(q)), potential function of q, or its gradient
#   Mf: Fisher observed(or expected) information matrix of approximating Normal
#   c: parameter to control step size of kappa
#   eps: step size
#   L: number of leapfrogs
# outputs:
#   q: new state of the parameter
#   u: new potential energy
#   Ind: proposal acceptance indicator

aMALA = function (q_cur, u_cur, U, Mf, c, eps=1)
{
  # initialization
  q = q_cur
  D = length(q)
  
  # sample kappa
  repeat
  {
    t=stats::runif(1,1/c,c)
    if(stats::runif(1)<(t+1/t)/(c+1/c))
      break
  }
  q[D]=q[D]*t
  
  # prepare pre-conditional matrix and gradient
  Q=Mf(q)
  cholQ=spam::chol(Q)
  g=U(q, grad = TRUE)$dlogpos
  
  # sample momentum
  z=stats::rnorm(D-1)
  p=spam::backsolve(cholQ,z)
  
  # log proposal density
  logprp = -t(z)%*%z/2+sum(log(spam::diag(cholQ)))
  
  # update momentum
  #	p=p-eps/2*solve(Q,g[-D])
  p = p-eps/2*(spam::chol2inv(cholQ)%*%g[-D])
  
  # update position
  q[-D] = q[-D]+eps*p
  
  # update pre-conditional matrix and gradient
  Q = Mf(c(q[-D],q_cur[D]))
  cholQ=spam::chol(Q)
  g = U(c(q[-D],q_cur[D]),T)$dlogpos # very interesting update!!!
  
  # update momentum
  p = p-eps/2*(spam::chol2inv(cholQ)%*%g[-D])
  
  # log reverse proposal density
  logprp_rev = -t(p)%*%Q%*%p/2+sum(log(spam::diag(cholQ)))
  
  # Evaluate potential energy
  pos_summ = U(q)
  u = pos_summ$logpos
  
  # Accept or reject the state jointly
  logAP = -u + u_cur - logprp + logprp_rev
  
  if ( is.finite(logAP) && (log(stats::runif(1))<min(0,logAP)) )
    return (list(q = q, u = u, Ind = 1, pos_summ = pos_summ))
  else
    return (list(q = q_cur, u = u_cur, Ind = 0, pos_summ = U(q_cur)))
}

#### Hamiltonian Monte Carlo ####
# This is standard HMC method.
# This function generates one sample given previous state

# inputs:
#   q_cur: initial state of the parameter
#   u_cur, du_cur: initial potential energy and its gradient
#   U:=-log(density(q)), potential function of q, or its gradient
#   eps: step size
#   L: number of leapfrogs
# outputs:
#   q: new state of the parameter
#   u, du: new potential energy and its gradient
#   Ind: proposal acceptance indicator

HMC = function (q_cur, u_cur, du_cur, U, eps=.2, L=5, rand_leap=TRUE)
{  
  # initialization
  q = q_cur
  D = length(q)
  u = u_cur
  du = du_cur
  
  # sample momentum
  p = stats::rnorm(D)
  
  # calculate current energy
  E_cur = u + sum(p^2)/2
  
  # Make a half step for momentum at the beginning
  p = p - eps/2 * du
  
  if (rand_leap)
    randL = ceiling(stats::runif(1)*L)
  else
    randL = ceiling(L)
  
  # Alternate full steps for position and momentum
  for (l in 1:randL)
  {
    # Make a full step for the position
    q = q + eps * p
    
    du = U(q, grad = TRUE)$dlogpos
    # Make a full step for the momentum, except at end of trajectory
    if (l!=randL)
      p = p - eps * du
  }
  
  # Make a half step for momentum at the end.
  p = p - eps/2 * du
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  pos_summ = U(q)
  u = pos_summ$logpos
  E_prp = u + sum(p^2)/2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  logAP = -E_prp + E_cur
  
  if(is.finite(logAP) && (log(stats::runif(1))<min(0,logAP)))
    return (list(q = q, u = u, du = du, Ind = 1, pos_summ = pos_summ))
  else
    return (list(q = q_cur, u = u_cur, du = du_cur, Ind = 0, pos_summ = U(q_cur)))
}

#### Split Hamiltonian Monte Carlo ####
# This is splitHMC method by (Gaussian) approximation.
# This function generates one sample given previous state.

# inputs:
#   q_cur: initial state of the parameter
#   u_cur, du_cur: initial potential energy and its gradient
#   U:=-log(density(q)), potential function of q, or its gradient
#   rtEV, EVC: square root of eigen-valudes, eigen-vectors of Fisher observed(or expected) information matrix of approximating Normal
#   eps: step size
#   L: number of leapfrogs
# outputs:
#   q: new state of the parameter
#   u, du: new potential energy and its gradient
#   Ind: proposal acceptance indicator

splitHMC = function (q_cur, u_cur, du_cur, U, rtEV, EVC, eps=.1, L=5, rand_leap=TRUE)
{
  # initialization
  q = q_cur
  D = length(q)
  u = u_cur
  du = du_cur
  
  # sample momentum
  p = stats::rnorm(D)
  
  # calculate current energy
  E_cur = u + sum(p^2)/2
  
  
  if (rand_leap)
    randL = ceiling(stats::runif(1)*L)
  else
    randL = ceiling(L)
  
  p = p - eps/2*du
  qT = rtEV*(t(EVC)%*%q[-D])
  pT = t(EVC)%*%p[-D]
  A = t(qT)%*%qT
  # Alternate full steps for position and momentum
  for (l in 1:randL)
  {
    p[D] <- p[D] - eps/2*A/2*exp(q[D])
    q[D] <- q[D] + eps/2*p[D]
    
    # Make a full step for the middle dynamics
    Cpx = complex(modulus = 1, argument = -rtEV*exp(q[D]/2)*eps)*complex(real = qT*exp(q[D]/2), imaginary = pT)
    qT = Re(Cpx)*exp(-q[D]/2)
    pT = Im(Cpx)
    q[-D] = EVC%*%(qT/rtEV)
    
    # Make a half step for the last half dynamics
    A=t(qT)%*%qT
    
    q[D] <- q[D] + eps/2*p[D]
    p[D] <- p[D] - eps/2*A/2*exp(q[D])
    
    du = U(q, grad = TRUE)
    if(l!=randL)
    {
      pT = pT - eps*(t(EVC)%*%du[-D])
      p[D] = p[D] - eps*du[D]
    }
  }
  p[-D] = EVC%*%pT - eps/2*du[-D]
  p[D] = p[D] - eps/2*du[D]
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  pos_summ = U(q)
  u = pos_summ$logpos
  E_prp = u + sum(p^2)/2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  logAP = -E_prp + E_cur
  
  if( is.finite(logAP) && (log(stats::runif(1))<min(0,logAP)) )
    return (list(q = q, u = u, du = du, Ind = 1, pos_summ = pos_summ))
  else
    return (list(q = q_cur, u = u_cur, du = du_cur, Ind = 0, pos_summ = U(q_cur)))
}

#### Helper functions ####

# Normal log-prior
log_mvnorm_prior <- function(x, prec, mu=rep(0, length(x)))
{
  x_mu = x - mu
  return(0.5 * spam::determinant(prec)$modulus - 0.5 * crossprod(x_mu, prec %*% x_mu))
}

# Gamma log-prior for kappa
log_kappa_prior <- function(kappa, alpha, beta)
{
  return(stats::dgamma(x = kappa, shape = alpha, rate = beta, log = TRUE))
}

# Gamma log-prior for tau
log_tau_prior <- function(tau, alpha, beta)
{
  return(stats::dgamma(x = exp(tau), shape = alpha+1, rate = beta, log = TRUE))
}

# Normal log-prior for betas
log_betas_prior <- function(betas, betas_prec = diag(1/100, 1/100))
{
  return(log_mvnorm_prior(x = betas, prec = betas_prec))
}

# Computes a summary of the posterior
compute_pos_summ = function(samp_alg, loglikf, f, kappa, invC, alpha, beta, lik_init,
                            betas = NULL, betas_prec = diag(1/100, 1/100), covar_vals = NULL,
                            loglik = NULL, logpri = NULL)
{
  if (is.null(loglik))
  {
    if (samp_alg %in% c("fixed", "MH"))
      loglik = loglikf(f, lik_init, betas, covar_vals)
    else if (samp_alg == "ESS")
      loglik = loglikf(c(f, betas), lik_init, length(f)+1, covar_vals)
    else
      loglik = loglikf(f, lik_init)
  }
  
  if (is.null(logpri))
  {
    logpri = log_mvnorm_prior(x = f, prec = invC * kappa) +
      log_kappa_prior(kappa = kappa, alpha = alpha, beta = beta)
    if (samp_alg %in% c("MH", "ESS"))
    {
      logbetapri = log_betas_prior(betas, betas_prec)
      logpri = logpri + logbetapri
    }
  }
  
  logpos = logpri + loglik
  
  return(list(loglik = loglik, logpri = logpri, logpos = logpos))
}

# Intrinsic precision matrix
Q_matrix <- function(input, s_noise = 0, signal = 1)
{
  n2 <- length(input)
  diff1 <- diff(input)
  diff1[diff1==0] <- s_noise #correction for dividing over 0
  diff <- (1/(signal*diff1))
  
  Q<-spam::spam(0,n2,n2)  
  if (n2>2)
  {
    Q[cbind(seq(1,n2),seq(1,n2))] <- c(diff[1], diff[1:(n2-2)] + diff[2:(n2-1)],
                                       diff[n2-1]) + (1/signal)*rep(s_noise, n2)
  }
  else
  {
    Q[cbind(seq(1,n2),seq(1,n2))] <- c(diff[1],diff[n2-1])+(1/signal)*rep(s_noise,n2)
  }
  Q[cbind(seq(1,n2-1),seq(2,n2))] <- -diff[1:(n2-1)]
  Q[cbind(seq(2,n2),seq(1,n2-1))] <- -diff[1:(n2-1)]
  
  return(Q)
}

# backwards compatibility (deprecate soon)
Q.matrix = function(...)
{
  return(Q_matrix(...))
}

burnin_subsample = function(res, burnin = 0, subsample = 1)
{
  nsamp = dim(res$samp)[1]
  indices = seq(burnin+1, nsamp, by=subsample)
  samp = res$samp[indices, ]
  pos_summ = res$pos_summ[indices, ]
  
  cleaned_result = list(samp=samp, pos_summ = pos_summ, grid = res$grid)
}

update_burnin_subsample = function(res, burnin = 0, subsample = 1)
{
  Ngrid = res$Ngrid
  cleaned_res = burnin_subsample(res = res, burnin = burnin, subsample = subsample)
  
  logfmat = cleaned_res$samp[,1:(Ngrid-1)]
  if (res$alg == "ESS" && res$samp_alg %in% c("MH", "ESS"))
  {
    params = cleaned_res$samp[,Ngrid:(Ngrid+2)]
  }
  else
  {
    params = matrix(cleaned_res$samp[,Ngrid])
  }
  estimates = calculate_estimates(logfmat = logfmat, params = params, grid = res$grid)
  
  res$cleaned_res = cleaned_res
  res$estimates = estimates
  
  res$med = estimates$fmed
  res$low = estimates$flow
  res$hi = estimates$fhi
  
  res$med_fun = estimates$fmed_fun
  res$low_fun = estimates$flow_fun
  res$hi_fun = estimates$fhi_fun
  
  return(res)
}

calculate_estimates = function(logfmat, params, grid)
{
  logfmed = apply(logfmat, MARGIN = 2, FUN = stats::median)
  logflow = apply(logfmat, MARGIN = 2, FUN = function(x) stats::quantile(x, .025))
  logfhi  = apply(logfmat, MARGIN = 2, FUN = function(x) stats::quantile(x, .975))
  
  logfmed_fun = stats::stepfun(grid, c(logfmed[1], logfmed, utils::tail(logfmed,1)))
  logflow_fun = stats::stepfun(grid, c(logflow[1], logflow, utils::tail(logflow,1)))
  logfhi_fun  = stats::stepfun(grid, c(logfhi[1], logfhi,  utils::tail(logfhi,1)))
  
  fmed = apply(exp(logfmat), MARGIN = 2, FUN = stats::median)
  flow = apply(exp(logfmat), MARGIN = 2, FUN = function(x) stats::quantile(x, .025))
  fhi  = apply(exp(logfmat), MARGIN = 2, FUN = function(x) stats::quantile(x, .975))
  
  fmed_fun = stats::stepfun(grid, c(fmed[1], fmed, utils::tail(fmed,1)))
  flow_fun = stats::stepfun(grid, c(flow[1], flow, utils::tail(flow,1)))
  fhi_fun  = stats::stepfun(grid, c(fhi[1], fhi,  utils::tail(fhi,1)))
  
  pmed = apply(params, MARGIN = 2, FUN = stats::median)
  plow = apply(params, MARGIN = 2, FUN = function(x) stats::quantile(x, .025))
  phi  = apply(params, MARGIN = 2, FUN = function(x) stats::quantile(x, .975))
  
  return(list(logfmed = logfmed, logflow = logflow, logfhi = logfhi,
              logfmed_fun = logfmed_fun, logflow_fun = logflow_fun,
              logfhi_fun = logfhi_fun, fmed = fmed, flow = flow, fhi = fhi,
              fmed_fun = fmed_fun, flow_fun = flow_fun, fhi_fun = fhi_fun,
              pmed = pmed, plow = plow, phi = phi))
}


#### Sampling functions ####

#MH_betas = function(curr_betas, lik_init, f, curr_pos_summ, proposal_sds = c(0.1, 0.1))
MH_betas = function(curr_betas, curr_pos_summ, loglikf, lik_init, f, kappa,
                    invC, alpha, beta, betas_prec = diag(1/100, 1/100),
                    covar_vals = covar_vals, proposal_sds = c(0.1, 0.1))
{
  # Metropolis step for betas
  new_betas = stats::rnorm(n = length(curr_betas), mean = curr_betas, sd = proposal_sds)
  
  new_pos_summ = compute_pos_summ(samp_alg = "MH", loglikf = loglikf, f = f,
                                  kappa = kappa, invC = invC, lik_init = lik_init,
                                  alpha = alpha, beta = beta, betas = new_betas,
                                  betas_prec = betas_prec, covar_vals = covar_vals)
  #new_ll = coal_samp_loglik(init = lik_init, f = f, beta0 = new_beta0, beta1 = new_beta1)
  
  if (new_pos_summ$logpos > curr_pos_summ$logpos ||
      log(stats::runif(n = 1)) < new_pos_summ$logpos - curr_pos_summ$logpos)
  {
    result = list(betas = new_betas, pos_summ = new_pos_summ, ind = 1)
  }
  else
  {
    result = list(betas = curr_betas, pos_summ = curr_pos_summ, ind = 0)
  }
  return(result)
}

MH_betas_rscan = function(curr_betas, curr_pos_summ, loglikf, lik_init, f, kappa,
                          invC, alpha, beta, betas_prec = diag(1/100, 1/100),
                          covar_vals = covar_vals, proposal_sds = rep(0.1, length(curr_betas)),
                          niter = length(curr_betas))
{
  n = length(curr_betas)
  rscan = sample(x = 1:n, size = niter, replace = TRUE)
  inds = vector(mode="list", length=n)
  
  new_betas = curr_betas
  for (idx in rscan)
  {
    new_betas[idx] = stats::rnorm(n = 1, mean = curr_betas[idx], sd = proposal_sds[idx])
    
    new_pos_summ = compute_pos_summ(samp_alg = "MH", loglikf = loglikf, f = f,
                                    kappa = kappa, invC = invC, lik_init = lik_init,
                                    alpha = alpha, beta = beta, betas = new_betas,
                                    betas_prec = betas_prec, covar_vals = covar_vals)
    
    if (new_pos_summ$logpos > curr_pos_summ$logpos ||
        log(stats::runif(n = 1)) < new_pos_summ$logpos - curr_pos_summ$logpos)
    {
      curr_betas = new_betas
      curr_pos_summ = new_pos_summ
      inds[[idx]] = c(inds[[idx]], 1)
    }
    else
    {
      new_betas = curr_betas
      inds[[idx]] = c(inds[[idx]], 0)
    }
  }
  
  result = list(betas = curr_betas, pos_summ = curr_pos_summ, inds = inds)
  
  return(result)
}

whiten_kappa = function(kappa, f, lik_init, cholC, invtcholC, loglikf, u, alpha, beta, prop_sd = 1)
{
  nu = (invtcholC * sqrt(kappa)) %*% f
  new_kappa = stats::rlnorm(n = 1, meanlog = log(kappa), sdlog = prop_sd)
  new_f = (t(cholC) / sqrt(new_kappa)) %*% nu
  
  new_u = coal_loglik(init = lik_init, f = new_f) # plus other terms
  udiff = new_u - u
  priordiff = log_kappa_prior(kappa = new_kappa, alpha = alpha, beta = beta) -
    log_kappa_prior(kappa = kappa, alpha = alpha, beta = beta)
  propdiff = stats::dlnorm(x = log(kappa), meanlog = log(new_kappa), sdlog = prop_sd) - 
    stats::dlnorm(x = log(new_kappa), meanlog = log(kappa), sdlog = prop_sd)
  
  if (log(stats::runif(n = 1)) < udiff + priordiff + propdiff)
  {
    result = list(kappa = new_kappa, f = new_f)
  }
  else
  {
    result = list(kappa = kappa, f = f)
  }
  
  return(result)
}

# This serves as a black box to sample distributions using HMC algorithms provided data and basic settings. #
sampling = function(data, para, alg, setting, init, verbose=TRUE, printevery=100)
{
  # pass the data and parameters
  lik_init = data$lik_init # f_offset = data$f_offset
  Ngrid = lik_init$ng+1
  alpha = para$alpha
  beta = para$beta
  invC = para$invC
  rtEV = para$rtEV
  EVC = para$EVC
  cholC = para$cholC
  
  # MCMC sampling setting
  stepsz = setting$stepsz
  Nleap  = setting$Nleap
  if (alg=='aMALA')
    szkappa = setting$szkappa
  
  if (alg=="HMC" | alg == "splitHMC")
  {
    rand_leap = setting$rand_leap
  }
  
  # if (alg == "ESSwS")
  # {
  #   betas = para$betas
  # }
  
  # storage of posterior samples
  NSAMP = setting$NSAMP
  NBURNIN = setting$NBURNIN
  NSUBSAMP = setting$NSUBSAMP
  recorded_iters = seq.int(from = NBURNIN+1, to = NSAMP, by = NSUBSAMP)
  
  samp = matrix(NA, length(recorded_iters), Ngrid) # all parameters together
  acpi = 0
  acpt = rep(NA, length(recorded_iters))
  
  # storage of log prior, log likelihood, and log posterior
  logpri = rep(NA, length(recorded_iters))
  loglik = rep(NA, length(recorded_iters))
  logpos = rep(NA, length(recorded_iters))
  
  # initialization
  theta = init$theta
  u = init$u
  du = init$du
  
  # if (alg %in% c("ESSwS-MH", "ESSwS-ESS"))
  # {
  #   betas = init$betas
  #   betas_out = NULL
  # }
  
  if (alg == "HMC")
  {
    Ufun = function(theta, grad=FALSE) U(theta = theta, init = lik_init, invC = invC,
                                         alpha = alpha, beta = beta, grad = grad)
    colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "tau")
  }
  else if (alg == "splitHMC")
  {
    Ufun = function(theta, grad=FALSE) U_split(theta = theta, init = lik_init, invC = invC,
                                               alpha = alpha, beta = beta, grad = grad)
    colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "tau")
  }
  else if (alg == "MALA")
  {
    Ufun = function(theta, grad=FALSE) U(theta = theta, init = lik_init, invC = invC,
                                         alpha = alpha, beta = beta, grad = grad)
    colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "tau")
  }
  else if (alg == "aMALA")
  {
    Ufun = function(theta, grad=FALSE) U_kappa(theta = theta, init = lik_init, invC = invC,
                                               alpha = alpha, beta = beta, grad = grad)
    Mf = function(theta) Met(theta, lik_init, invC)
    colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "kappa")
  }

  # start MCMC run
  start_time = Sys.time()
  cat('Running ', alg ,' sampling...\n')
  for(iter in 1:NSAMP)
  {
    if (alg == "HMC")
    {
      #Ufun = function(theta, grad=FALSE) U(theta = theta, init = lik_init, invC = invC,
      #                                     alpha = alpha, beta = beta, grad = grad)
      res = HMC(theta, u, du, Ufun, stepsz, Nleap, rand_leap)
    }
    else if (alg == "splitHMC")
    {
      #Ufun = function(theta, grad=FALSE) U_split(theta = theta, init = lik_init, invC = invC,
      #                                           alpha = alpha, beta = beta, grad = grad)
      res = splitHMC(theta, u, du, Ufun, rtEV, EVC, stepsz, Nleap, rand_leap)
    }
    else if (alg == "MALA")
    {
      #Ufun = function(theta, grad=FALSE) U(theta = theta, init = lik_init, invC = invC,
      #                                     alpha = alpha, beta = beta, grad = grad)
      res = MALA(theta, u, du, Ufun, stepsz)
    }
    else if (alg == "aMALA")
    {
      #Ufun = function(theta, grad=FALSE) U_kappa(theta = theta, init = lik_init, invC = invC,
      #                                           alpha = alpha, beta = beta, grad = grad)
      #Mf = function(theta) Met(theta, lik_init, invC)
      res = aMALA(q_cur = theta, u_cur = u, U = Ufun, Mf = Mf, c = szkappa, eps = stepsz)
    }
    else
    {
      stop('The algorithm is not in the list!')
    }
    
    acpi <- acpi+res$Ind
    
    theta[1:Ngrid] <- res$q

    u <- res$u
    if (alg %in% c('HMC','splitHMC','MALA'))
      du <- res$du
    
    # save posterior samples after burnin
    output_index = match(x = iter, table = recorded_iters)
    if (!is.na(output_index))
    {
      samp[output_index, ] <- theta
      acpt[output_index] <- res$Ind
      
      logpri[output_index] <- res$pos_summ$logpri
      loglik[output_index] <- res$pos_summ$loglik
      logpos[output_index] <- res$pos_summ$logpos
    }
    
    if(verbose && iter %% printevery == 0)
    {
      cat(iter, ' iterations have been finished!\n' )
      cat('Online acceptance rate is ',acpi/printevery,'\n')
      acpi=0
    }
  }
  stop_time <- Sys.time()
  time <- stop_time-start_time
  cat('\nTime consumed : ', time, units(time))
  #acpt <- acpt/(NSAMP-NBURNIN)
  cat('\nFinal Acceptance Rate: ', sum(acpt) / (NSAMP-NBURNIN),'\n')
  
  pos_summ = data.frame(acpt=acpt, logpri = logpri, loglik = loglik, logpos = logpos)
  
  result = list(samp=samp, alg=alg, time=time, pos_summ = pos_summ)

  return(result)
}

ESS_none_ll = function(f, lik_init)
{
  return(coal_loglik(init = lik_init, f = f))
}

ESS_betas_ll = function(f, lik_init, betas, covar_vals = NULL) 
{
  fs = cbind(f, covar_vals, deparse.level = 0)
  return(coal_loglik(init = lik_init, f = f) +
           samp_loglik(init = lik_init, fs = fs, betas = betas))
}

ESS_ext_ll = function(f, lik_init, Ngrid, covar_vals = NULL) 
{
  nbetas = length(f) - Ngrid + 1
  return(coal_loglik(init = lik_init, f = f[1:(Ngrid-1)]) +
           samp_loglik(init = lik_init,
                       fs = cbind(f[1:(Ngrid-1)], covar_vals, deparse.level = 0),
                       betas = f[Ngrid:(Ngrid+nbetas-1)]))
}

sampling_ESS = function(data, para, setting, init,
                        samp_alg = "none", kappa_alg = "gibbs",
                        verbose=TRUE, printevery=100)
{
  # pass the data and parameters
  lik_init = data$lik_init
  covar_vals = data$covar_vals
  Ngrid = lik_init$ng+1
  
  alpha = para$alpha
  beta = para$beta
  invC = para$invC
  cholC = para$cholC
  beta_vars = para$beta_vars
  
  proposal_sds = setting$proposal_sds
  
  # storage of posterior samples
  NSAMP = setting$NSAMP
  NBURNIN = setting$NBURNIN
  NSUBSAMP = setting$NSUBSAMP
  recorded_iters = seq.int(from = NBURNIN+1, to = NSAMP, by = NSUBSAMP)
  
  # initialization
  f = init$theta[1:(Ngrid-1)]
  kappa = init$theta[Ngrid]
  #u = init$u
  
  #acpi = 0
  acpt = rep(1, length(recorded_iters))
  
  fmat = matrix(NA, nrow = length(recorded_iters), ncol = length(f))
  kappas = rep(NA, length(recorded_iters))
  
  # storage of log prior, log likelihood, and log posterior
  logpri = rep(NA, length(recorded_iters))
  loglik = rep(NA, length(recorded_iters))
  logpos = rep(NA, length(recorded_iters))
  
  if (samp_alg == "none")
  {
    ll = ESS_none_ll
  }
  else if (samp_alg == "fixed")
  {
    betas = para$betas
    ll = ESS_betas_ll
  }
  else if (samp_alg == "MH")
  {
    betas = init$betas
    ll = ESS_betas_ll
  }
  else if (samp_alg == "ESS")
  {
    betas = init$betas
    ll = ESS_ext_ll
  }
  
  if (samp_alg %in% c("MH", "ESS"))
  {
    betas_out = matrix(NA, nrow = length(recorded_iters), ncol = length(betas))
  }
  
  pos_summ = compute_pos_summ(samp_alg = samp_alg, loglikf = ll, f = f,
                              kappa = kappa, invC = invC, lik_init = lik_init,
                              alpha = alpha, beta = beta, betas = betas,
                              betas_prec = diag(1/beta_vars), covar_vals = covar_vals)
  
  if (kappa_alg == "whiten")
  {
    invtcholC = solve(t(cholC))
  }
  
  # start MCMC run
  start_time = Sys.time()
  cat('Running ESS', samp_alg ,' sampling...\n')
  for(iter in 1:NSAMP)
  {
    if (samp_alg == "none")
    {
      res = ESS(q_cur = f, l_cur = pos_summ$loglik, loglik = ll,
                cholC = cholC/sqrt(kappa), lik_init = lik_init)
      f = res$q
    }
    else if (samp_alg == "fixed")
    {
      res = ESS(q_cur = f, l_cur = pos_summ$loglik, loglik = ll,
                cholC = cholC/sqrt(kappa), lik_init = lik_init,
                betas = betas, covar_vals = covar_vals)
      f = res$q
    }
    else if (samp_alg == "MH")
    {
      MH_res = MH_betas_rscan(curr_betas = betas, curr_pos_summ = pos_summ,
                              loglikf = ll, lik_init = lik_init, f = f, kappa = kappa,
                              invC = invC, alpha = alpha, beta = beta,
                              betas_prec = diag(1/beta_vars), covar_vals = covar_vals,
                              proposal_sds = proposal_sds)
      betas = MH_res$betas
      pos_summ = MH_res$pos_summ
      
      res = ESS(q_cur = f, l_cur = pos_summ$loglik, loglik = ll,
                cholC = cholC/sqrt(kappa), lik_init = lik_init,
                betas = betas, covar_vals = covar_vals)
      f = res$q
    }
    else if (samp_alg == "ESS")
    {
      cholC_betas = matrix(0, nrow = nrow(cholC)+2, ncol = ncol(cholC)+2)
      cholC_betas[1:nrow(cholC), 1:ncol(cholC)] = cholC/sqrt(kappa)
      cholC_betas[nrow(cholC)+1, ncol(cholC)+1] = sqrt(beta_vars[1])
      cholC_betas[nrow(cholC)+2, ncol(cholC)+2] = sqrt(beta_vars[2])
      
      res = ESS(q_cur = c(f, betas), l_cur = pos_summ$loglik, loglik = ll,
                cholC = cholC_betas, lik_init = lik_init,
                Ngrid = Ngrid, covar_vals = covar_vals)
      f = res$q[1:(Ngrid-1)]
      betas = res$q[Ngrid:(Ngrid+1)]
    }
    else
    {
      stop('The ESS subalgorithm is not in the list!')
    }
    
    #acpi <- acpi+res$Ind
    
    if (kappa_alg == "gibbs")
    {
      kappa <- stats::rgamma(1, alpha + (Ngrid-1)/2, beta + crossprod(f, invC %*% f)/2)
    }
    else if (kappa_alg == "whiten")
    {
      kappa_res <- whiten_kappa(kappa = kappa, f = f, lik_init = lik_init,
                                cholC = cholC, invtcholC = invtcholC,
                                loglikf = ll, u = pos_summ$loglik, alpha = alpha, beta = beta)
      kappa <- kappa_res$kappa
      f <- kappa_res$f
    }
    else
      stop("Kappa operator not recognized.")
    
    pos_summ = compute_pos_summ(samp_alg = samp_alg, loglikf = ll, f = f,
                                lik_init = lik_init, kappa = kappa, invC = invC,
                                alpha = alpha, beta = beta, betas = betas,
                                betas_prec = diag(1/beta_vars), covar_vals = covar_vals)
    #u <- pos_summ$loglik
    
    # save posterior samples after burnin
    output_index = match(x = iter, table = recorded_iters)
    if (!is.na(output_index))
    {
      fmat[output_index, ] <- f
      kappas[output_index] <- kappa
      
      if (samp_alg %in% c("MH", "ESS"))
      {
        betas_out[output_index, ] <- betas
      }
      #acpt[output_index] <- res$Ind
      
      logpri[output_index] <- pos_summ$logpri
      loglik[output_index] <- pos_summ$loglik
      logpos[output_index] <- pos_summ$logpos
    }
    
    if (verbose && iter %% printevery == 0)
    {
      cat(iter, ' iterations have been finished!\n' )
      cat('Online acceptance rate is ', 1,'\n')
      #acpi=0
    }
  }
  stop_time <- Sys.time()
  time <- stop_time-start_time
  cat('\nTime consumed : ', time, units(time))
  #acpt <- acpt/(NSAMP-NBURNIN)
  cat('\nFinal Acceptance Rate: ', sum(acpt) / length(recorded_iters),'\n')
  
  pos_summ = data.frame(acpt=acpt, logpri = logpri, loglik = loglik, logpos = logpos)
  
  if (samp_alg %in% c("MH", "ESS"))
  {
    samp = cbind(fmat, kappas, betas_out)
    colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "kappa", paste("beta", 0:(length(betas)-1)))
  }
  else
  {
    samp = cbind(fmat, kappas)
    colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "kappa")
  }
  
  return(list(samp=samp, alg="ESS", time=time, pos_summ = pos_summ,
              samp_alg = samp_alg, kappa_alg = kappa_alg))
}

#' MCMC Sampling
#' 
#' @param dataset \code{phylo} object or list containing vectors of coalescent 
#'   times \code{coal_times}, sampling times \code{samp_times}, and number 
#'   sampled per sampling time \code{n_sampled}.
#' @param alg string selecting which MCMC sampler to use. Options are "HMC", 
#'   "splitHMC", "MALA", "aMALA", and "ESS".
#' @param nsamp integer number of MCMC steps to compute.
#' @param nburnin integer number of MCMC steps to discard as burn-in.
#' @param nsubsamp integer after burn-in, how often to record a step to the 
#'   output.
#' @param ngrid integer number of grid point in the latent field.
#' @param nugget string selecting which "nugget" adjustment to apply to the 
#'   precision matrix to make it full-rank. Options are '1,1' for an adjustment 
#'   to the first element, 'diag' for an adjustment to the entire main diagonal,
#'   or 'none' which may result in a non-full-rank precision matrix.
#' @param prec_alpha,prec_beta numeric shape and rate parameters for the prior 
#'   on precision.
#' @param TrjL numeric tuning parameter.
#' @param Nleap integer tuning parameter.
#' @param szkappa numeric tuning parameter.
#' @param rand_leap logical tuning parameter.
#' @param f_init numeric vector starting log effective population size values.
#' @param kappa numeric starting kappa.
#' @param covariates list of functions representing covariate trajectories that 
#'   (may) influence sampling frequency.
#' @param betas numeric vector of starting values for the beta hyperparameters.
#' @param samp_alg string selecting sampling algorithm for sampling time 
#'   intensity coefficients. One of "none" (default), "fixed", "MH", and "ESS".
#' @param kappa_alg selects sampling algorithm for kappa. One of "gibbs" 
#'   (default) or "whiten".
#' @param beta_vars numeric vector prior variances of the beta hyperparameters.
#' @param printevery integer how many MCMC steps between writing output to the 
#'   console.
#'   
#' @export
mcmc_sampling = function(dataset, alg, nsamp, nburnin=0, nsubsamp=1, ngrid=100,
                         nugget="1,1", prec_alpha = 1e-2, prec_beta = 1e-2,
                         TrjL=NULL, Nleap=NULL, szkappa=NULL, rand_leap=NULL,
                         f_init = rep(1, ngrid-1), kappa = 1,
                         covariates=NULL, betas=rep(0, 2+length(covariates)),
                         samp_alg = "none", kappa_alg = "gibbs",
                         beta_vars = rep(100, length(betas)), printevery=100)
{
  if (class(dataset) == "phylo")
  {
    phy <- summarize_phylo(dataset)
  }
  else if (all(c("coal_times", "samp_times", "n_sampled") %in% names(dataset)))
  {
    phy <- with(dataset, list(samp_times = samp_times, coal_times = coal_times,
                           n_sampled = n_sampled))
  }
  
  samp_times = phy$samp_times
  n_sampled  = phy$n_sampled
  coal_times = phy$coal_times
  
  # Jump tuning parameters--should probably have an option to change in the arguments
  if (is.null(TrjL))
    TrjL = switch(alg, HMC=3, splitHMC=3, MALA=0.1, aMALA=0.1)
  if (is.null(Nleap))
    Nleap = switch(alg, HMC=30, splitHMC=15, MALA=1, aMALA=1)
  if (is.null(szkappa) & alg=="aMALA")
    szkappa = 1.2
  
  if (is.null(rand_leap) & (alg=="HMC" | alg=="splitHMC"))
    rand_leap = TRUE
  
  stepsz = TrjL/Nleap
  
  grid_bds = range(c(coal_times,samp_times))
  #Ngrid = 100
  
  grid = seq(grid_bds[1],grid_bds[2],length.out=ngrid)
  intl = grid[2]-grid[1]
  midpts = grid[-1]-intl/2
  
  covar_vals = NULL
  if (!is.null(covariates))
  {
    for (fcn in covariates)
    {
      covar_vals = cbind(covar_vals, log(fcn(midpts)), deparse.level = 0)
    }
  }
  
  # initialize likelihood calculation
  lik_init = coal_lik_init(samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times, grid=grid)
  
  # calculate intrinsic precision matrix
  invC <- Q_matrix(midpts,0,1)
  
  # fudge to be able to compute the cholC
  if (nugget == "1,1")
    invC[1,1] <- invC[1,1]+.0001 # nugget at (1,1)
  else if (nugget == "diag")
  #Julia: I had a warning when using this--needs to be corrected
    diag(invC)<-diag(invC)+.0001 # nugget for the whole diagonal
  else if (nugget == "none")
    warning("No nugget may result in a non-full-rank matrix.")
  else
    stop(paste("Unrecognized argument nugget = '", nugget, "', please use '1,1', 'diag', or 'none'.", sep = ""))
  
  eig  = spam::eigen.spam(invC, TRUE)
  rtEV = sqrt(eig$values)
  EVC  = eig$vectors
  
  C = spam::solve.spam(invC)
  cholC = chol(C)
  
  # initializations
  #theta = rep(1,Ngrid)
  theta = c(f_init, kappa)
  
  if (alg == "HMC")
  {
    u  = U(theta,lik_init,invC,prec_alpha,prec_beta)$logpos
    du = U(theta,lik_init,invC,prec_alpha,prec_beta, TRUE)$dlogpos
  }
  else if (alg == "splitHMC")
  {
    u  = U_split(theta,lik_init,invC,prec_alpha,prec_beta)$logpos
    du = U_split(theta,lik_init,invC,prec_alpha,prec_beta, TRUE)
  }
  else if (alg == "MALA")
  {
    u  = U(theta,lik_init,invC,prec_alpha,prec_beta)$logpos
    du = U(theta,lik_init,invC,prec_alpha,prec_beta, TRUE)$dlogpos
  }
  else if (alg == "aMALA")
  {
    u  = U_kappa(theta,lik_init,invC,prec_alpha,prec_beta)$logpos
    du = NULL
  }
  else if (alg == "ESS")
  {
    u = NULL
    #u  = coal_loglik(init = lik_init, f = theta[-Ngrid])
    du = NULL
  }
  else
  {
    stop('The algorithm is not in the list!')
  }
  
  # MCMC sampling preparation
  dataset = list(lik_init = lik_init, covar_vals = covar_vals)
  para = list(alpha = prec_alpha, beta = prec_beta, invC = invC, rtEV = rtEV,
              EVC = EVC, cholC = cholC, betas = betas, beta_vars = beta_vars)
  setting = list(stepsz = stepsz, Nleap = Nleap,
                 NSAMP = nsamp, NBURNIN = nburnin, NSUBSAMP = nsubsamp,
                 szkappa = szkappa, rand_leap=rand_leap,
                 proposal_sds = rep(0.3, length(betas)))
  init = list(theta = theta, u = u, du = du, betas = betas)
  
  # Run MCMC sampler
  if (alg == "ESS")
  {
    res_MCMC = sampling_ESS(data = dataset, para = para, setting = setting,
                            init = init, samp_alg = samp_alg, kappa_alg = kappa_alg,
                            printevery = printevery)
  }
  else
  {
    res_MCMC = sampling(data = dataset, para = para, alg = alg, setting = setting,
                        init = init, printevery = printevery)
  }
  
  res_MCMC$alg = alg
  res_MCMC$samp_alg = samp_alg
  res_MCMC$kappa_alg = kappa_alg
  res_MCMC$Ngrid = ngrid
  
  #cleaned_res = burnin_subsample(res = res_MCMC, burnin = 0)
  
  logfmat = res_MCMC$samp[,1:(ngrid-1)]
  if (alg == "ESS" && samp_alg %in% c("MH", "ESS"))
  {
    params = res_MCMC$samp[,ngrid:(ngrid+2)]
  }
  else
  {
    params = matrix(res_MCMC$samp[,ngrid])
  }
  estimates = calculate_estimates(logfmat = logfmat, params = params, grid = grid)
  
  #res_MCMC$cleaned_res = cleaned_res
  res_MCMC$estimates = estimates
  
  res_MCMC$med = estimates$fmed
  res_MCMC$low = estimates$flow
  res_MCMC$hi = estimates$fhi
  
  res_MCMC$med_fun = estimates$fmed_fun
  res_MCMC$low_fun = estimates$flow_fun
  res_MCMC$hi_fun = estimates$fhi_fun
  
  res_MCMC$grid = grid
  res_MCMC$x = midpts
  res_MCMC$samp_times = samp_times
  res_MCMC$n_sampled = n_sampled
  res_MCMC$coal_times = coal_times
  
  return(res_MCMC)
}

#' SMC' sampler--inference from local genealogies
#' 
#' @param data a list containing sufficient statistics
#' @param nsamp integer specifying number of MCMC samples
#' @param nburnin integer specifying the number of burnin samples
#' @param grid a vector with the grid points
#' @param alpha hyperparameter of precision of BM prior
#' @param beta hyperparameter of precision of BM prior
#' @param stepsz numeric tuning parameter for Split Hamiltonian Monte Carlo
#' @param Nleap integer tuning parameter for Split Hamiltonian Monte Carlo
#' @param rand_leap tuning parameter for Split Hamiltonian Monte Carlo
#' @param scaling numeric re-scaling parameter
#' @param tol tolerance to detect difference
#' 
#' @return A matrix of sim rows. Entry x_{i,j} has the n-j+1-th coalescent time of the i-th tree
#' @export
smcp_sampling = function(data,  nsamp, nburnin, grid, alpha = 1e-3, beta = 1e-3,
                         stepsz=.1, Nleap=15, rand_leap=TRUE,scaling=10,tol=1e-5)
{
  Ngrid<-length(grid)-1
  
  #MCMC sampling preparation
  SAMP = list(2)
  SAMP[[1]] = matrix(NA,nsamp-nburnin,Ngrid) # transformed effective population size
  SAMP[[2]] = rep(NA,nsamp-nburnin) # precision parameter in Brownian motion
  acpi = 0
  acpt = 0
  PRINT = TRUE
  
  #Initial Values
  f_init = rep(0.5,Ngrid)
  theta <- c(log(f_init),-1.6)+.0001
  alldata <- get_data(grid,data$sim,data$D,data$n,data$info_times,data$Fl,data$latent,data$t_new,data$t_del,tol)
  
  U<-function(theta,grad=F)U_split_smc(theta,alldata$lik_init,alldata$invC,alpha,beta,grad)
  current.u<-U(theta,F)$logpos
  current.grad<-U(theta,T)
  
  for(Iter in 1:nsamp)
  {
    if(PRINT&&Iter%%50==0)
    {
      cat(Iter, ' iterations have been finished!\n' )
      cat('Online acceptance rate is ',acpi/50,'\n')
      acpi=0
    }
    res=eval(parse(text='splitHMC'))(theta,current.u,current.grad,function(theta,grad=F)U_split_smc(theta,alldata$lik_init,alldata$invC,alpha,beta,grad),alldata$rtEV,alldata$EVC,stepsz,Nleap,rand_leap)
    theta=res$q;
    current.u<-res$u
    current.grad<-res$du
    N<-exp(theta[1:(length(theta)-1)])
    acpi=acpi+res$Ind
    if(Iter>nburnin)
    {
      SAMP[[1]][Iter-nburnin,]<-theta[1:(length(theta)-1)]
      SAMP[[2]][Iter-nburnin]<-theta[length(theta)]
      acpt<-acpt+res$Ind
    }
  }
  ini<-1
  med=apply(SAMP[[1]][ini:(Iter-nburnin-1),], MARGIN = 2, FUN = stats::median);
  low=apply(SAMP[[1]][ini:(Iter-nburnin-1),], MARGIN = 2, FUN = function(x)stats::quantile(x,.025))
  up=apply(SAMP[[1]][ini:(Iter-nburnin-1),], MARGIN = 2, FUN = function(x)stats::quantile(x,.975))
  
  results<-cbind(grid/scaling,c(low[1]-log(scaling),low-log(scaling)),c(med[1]-log(scaling),med-log(scaling)),c(up[1],up)-log(scaling))
  return(results)
}


