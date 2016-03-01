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

ESS = function(q_cur, l_cur, loglik, cholC)
{  
  # choose ellipse
  nu = t(cholC) %*% rnorm(length(q_cur))
  
  # log-likelihood threshold
  u = runif(1)
  logy <- l_cur + log(u)
  
  # draw a initial proposal, also defining a bracket
  t = 2*pi*runif(1)
  t_min <- t-2*pi
  t_max <- t
  
  q <- q_cur*cos(t) + nu*sin(t)
  l <- loglik(q)
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
    
    t <- runif(1, t_min, t_max)
    q <- q_cur*cos(t) + nu*sin(t)
    l <- loglik(q)
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
  nu = t(cholC) %*% rnorm(length(q_cur))
  
  # log-likelihood threshold
  u = runif(1)
  logy <- l_cur + log(u)
  
  # draw a initial proposal, also defining a bracket
  t = 2*pi*runif(1)
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
    
    t = runif(1, t_min, t_max)
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
  p = rnorm(D)
  
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
  
  if( is.finite(logAP) && (log(runif(1)) < min(0, logAP)) ) 
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
    t=runif(1,1/c,c)
    if(runif(1)<(t+1/t)/(c+1/c))
      break
  }
  q[D]=q[D]*t
  
  # prepare pre-conditional matrix and gradient
  Q=Mf(q)
  cholQ=spam::chol.spam(Q)
  g=U(q, grad = TRUE)$dlogpos
  
  # sample momentum
  z=rnorm(D-1)
  p=spam::backsolve.spam(cholQ,z)
  
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
  
  if ( is.finite(logAP) && (log(runif(1))<min(0,logAP)) )
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
  p = rnorm(D)
  
  # calculate current energy
  E_cur = u + sum(p^2)/2
  
  # Make a half step for momentum at the beginning
  p = p - eps/2 * du
  
  if (rand_leap)
    randL = ceiling(runif(1)*L)
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
  
  if(is.finite(logAP) && (log(runif(1))<min(0,logAP)))
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
  p = rnorm(D)
  
  # calculate current energy
  E_cur = u + sum(p^2)/2
  
  
  if (rand_leap)
    randL = ceiling(runif(1)*L)
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
    Cpx = complex(mod=1,arg=-rtEV*exp(q[D]/2)*eps)*complex(re=qT*exp(q[D]/2),im=pT)
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
  
  if( is.finite(logAP) && (log(runif(1))<min(0,logAP)) )
    return (list(q = q, u = u, du = du, Ind = 1, pos_summ = pos_summ))
  else
    return (list(q = q_cur, u = u_cur, du = du_cur, Ind = 0, pos_summ = U(q_cur)))
}

#### Helper functions ####

# Normal log-prior
log_mvnorm_prior <- function(x, prec, mu=rep(0, length(x)))
{
  return(-0.5 * t(x-mu) %*% prec %*% (x-mu))
}

# Gamma log-prior for kappa
log_kappa_prior <- function(kappa, alpha, beta)
{
  return(dgamma(x = kappa, shape = alpha, rate = beta, log = TRUE))
}

# Gamma log-prior for tau
log_tau_prior <- function(tau, alpha, beta)
{
  return(dgamma(x = exp(tau), shape = alpha+1, rate = beta, log = TRUE))
}

# Normal log-prior for betas
log_betas_prior <- function(betas, betas_prec = diag(1/100, 1/100))
{
  return(log_mvnorm_prior(x = betas, prec = betas_prec))
}

# Intrinsic precision matrix
Q_matrix <- function(input, s_noise, signal)
{
  n2 <- nrow(input)
  diff1 <- diff(input)
  diff1[diff1==0] <- s_noise #correction for dividing over 0
  diff <- (1/(signal*diff1))
  
  Q<-spam::spam(0,n2,n2)  
  if (n2>2)
  {
    Q[cbind(seq(1,n2),seq(1,n2))] <- c(diff[1], diff[1:(n2-2)] + diff[2:(n2-1)], diff[n2-1]) + (1/signal)*rep(s_noise, n2)
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

#### Sampling wrappers ####

# This serves as a black box to sample distributions using HMC algorithms provided data and basic settings. #
sampling = function(data, para, alg, setting, init, verbose=TRUE)
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
  
  if (alg == "ESSwS")
  {
    betas = para$betas
  }
  
  # storage of posterior samples
  NSAMP = setting$NSAMP
  #NBURNIN = setting$NBURNIN
  NBURNIN = 0
  
  samp = matrix(NA, NSAMP - NBURNIN, Ngrid) # all parameters together
  acpi = 0
  acpt = rep(NA, NSAMP - NBURNIN)
  
  # storage of log prior, log likelihood, and log posterior
  logpri = rep(NA, NSAMP - NBURNIN)
  loglik = rep(NA, NSAMP - NBURNIN)
  logpos = rep(NA, NSAMP - NBURNIN)
  
  # initialization
  theta = init$theta
  u = init$u
  du = init$du
  
  if (alg %in% c("ESSwS-MH", "ESSwS-ESS"))
  {
    betas = init$betas
    betas_out = NULL
  }
  
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
  else if (alg == "ESS")
  {
    ll = function(f) coal_loglik(init = lik_init, f = f)
    colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "kappa")
  }
  else if (alg == "ESSwS")
  {
    ll = function(f) coal_samp_loglik(init = lik_init, f = f, beta0 = betas[1], beta1 = betas[2])
    colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "kappa")
  }
  else if (alg == "ESSwS-MH")
  {
    ll = function(f) coal_samp_loglik(init = lik_init, f = f, beta0 = betas[1], beta1 = betas[2])
    colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "kappa")
  }
  else if (alg == "ESSwS-ESS")
  {
    ll = function(f) coal_samp_loglik(init = lik_init, f = f[1:(Ngrid-1)], beta0 = f[Ngrid], beta1 = f[Ngrid+1])
    colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "kappa")
  }
  
  # start MCMC run
  start_time = Sys.time()
  cat('Running ', alg ,' sampling...\n')
  for(Iter in 1:NSAMP)
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
    else if (alg == "ESS")
    {
      #ll = function(f) coal_loglik(init = lik_init, f = f)
      res = ESS(q_cur = theta[-Ngrid], l_cur = u, loglik = ll, cholC = cholC/sqrt(theta[Ngrid]))
      pos_summ = list(loglik = res$u)
      pos_summ$logpri = log_mvnorm_prior(x = res$q, prec = invC * theta[Ngrid]) +
                       log_kappa_prior(kappa = theta[Ngrid], alpha = alpha, beta = beta)
      pos_summ$logpos = pos_summ$logpri + pos_summ$loglik
      res$pos_summ = pos_summ
    }
    else if (alg == "ESSwS")
    {
      #ll = function(f) coal_samp_loglik(init = lik_init, f = f, beta0 = betas[1], beta1 = betas[2])
      res = ESS(q_cur = theta[-Ngrid], l_cur = u, loglik = ll, cholC = cholC/sqrt(theta[Ngrid]))
      pos_summ = list(loglik = res$u)
      pos_summ$logpri = log_mvnorm_prior(x = res$q, prec = invC * theta[Ngrid]) +
                       log_kappa_prior(kappa = theta[Ngrid], alpha = alpha, beta = beta)
      pos_summ$logpos = pos_summ$logpri + pos_summ$loglik
      res$pos_summ = pos_summ
    }
    else if (alg == "ESSwS-MH")
    {
      # Metropolis step for betas
      new_beta0 = rnorm(n = 1, mean = betas[1], sd = 0.1)
      new_beta1 = rnorm(n = 1, mean = betas[2], sd = 0.1)
      new_u = coal_samp_loglik(init = lik_init, f = theta[-Ngrid], beta0 = new_beta0, beta1 = new_beta1)
      if (new_u > u || log(runif(n = 1)) < new_u - u)
      {
        u = new_u
        betas = c(new_beta0, new_beta1)
      }
      
      #ll = function(f) coal_samp_loglik(init = lik_init, f = f, beta0 = betas[1], beta1 = betas[2])
      res = ESS(q_cur = theta[-Ngrid], l_cur = u, loglik = ll, cholC = cholC/sqrt(theta[Ngrid]))
      pos_summ = list(loglik = res$u)
      pos_summ$logpri = log_mvnorm_prior(x = res$q, prec = invC * theta[Ngrid]) +
                       log_kappa_prior(kappa = theta[Ngrid], alpha = alpha, beta = beta)
      pos_summ$logpos = pos_summ$logpri + pos_summ$loglik
      res$pos_summ = pos_summ
    }
    else if (alg == "ESSwS-ESS")
    {
      cholC_betas = matrix(0, nrow = nrow(cholC)+2, ncol = ncol(cholC)+2)
      cholC_betas[1:nrow(cholC), 1:ncol(cholC)] = cholC/sqrt(theta[Ngrid])
      cholC_betas[nrow(cholC)+1, ncol(cholC)+1] = sqrt(100)
      cholC_betas[nrow(cholC)+2, ncol(cholC)+2] = sqrt(100)
      
      #ll = function(f) coal_samp_loglik(init = lik_init, f = f[1:(Ngrid-1)], beta0 = f[Ngrid], beta1 = f[Ngrid+1])
      res = ESS(q_cur = c(theta[-Ngrid], betas), l_cur = u, loglik = ll, cholC = cholC_betas)
      pos_summ = list(loglik = res$u)
      pos_summ$logpri = log_mvnorm_prior(x = res$q[1:(Ngrid-1)], prec = invC * theta[Ngrid]) +
        log_mvnorm_prior(x = res$q[Ngrid:(Ngrid+1)], prec = diag(c(1/100, 1/100))) +
        log_kappa_prior(kappa = theta[Ngrid], alpha = alpha, beta = beta)
      pos_summ$logpos = pos_summ$logpri + pos_summ$loglik
      res$pos_summ = pos_summ
    }
    else
    {
      stop('The algorithm is not in the list!')
    }
    
    acpi <- acpi+res$Ind
    
    if (alg %in% c("ESS", "ESSwS", "ESSwS-MH"))
    {
      theta[1:(Ngrid-1)] <- res$q
      
      # Gibbs sample kappa for ESS
      theta[Ngrid] <- rgamma(1,alpha+(Ngrid-1)/2, beta + t(theta[-Ngrid]) %*% invC %*% theta[-Ngrid]/2)
    }
    else if (alg == "ESSwS-ESS")
    {
      theta[1:(Ngrid-1)] <- res$q[1:(Ngrid-1)]
      betas = res$q[Ngrid:(Ngrid+1)]
      
      # Gibbs sample kappa for ESS
      theta[Ngrid] <- rgamma(1,alpha+(Ngrid-1)/2, beta + t(theta[1:(Ngrid-1)]) %*% invC %*% theta[1:(Ngrid-1)]/2)
    }
    else
    {
      theta[1:Ngrid] <- res$q
    }
    #theta[1:(Ngrid-(alg=='ESS'))]=res$q
    
    u <- res$u
    if (alg %in% c('HMC','splitHMC','MALA'))
      du <- res$du
    
    # save posterior samples after burnin
    if (Iter > NBURNIN)
    {
      samp[Iter - NBURNIN, ] <- theta
      acpt[Iter - NBURNIN] <- res$Ind
      
      logpri[Iter - NBURNIN] <- res$pos_summ$logpri
      loglik[Iter - NBURNIN] <- res$pos_summ$loglik
      logpos[Iter - NBURNIN] <- res$pos_summ$logpos
      
      if (alg %in% c("ESSwS-MH", "ESSwS-ESS"))
      {
        betas_out = rbind(betas_out, betas, deparse.level = 0)
      }
    }
    
    if(verbose && Iter %% 100 == 0)
    {
      cat(Iter, ' iterations have been finished!\n' )
      cat('Online acceptance rate is ',acpi/100,'\n')
      acpi=0
    }
  }
  stop_time <- Sys.time()
  time <- stop_time-start_time
  cat('\nTime consumed : ',time)
  #acpt <- acpt/(NSAMP-NBURNIN)
  cat('\nFinal Acceptance Rate: ', sum(acpt) / (NSAMP-NBURNIN),'\n')
  
  pos_summ = data.frame(acpt=acpt, logpri = logpri, loglik = loglik, logpos = logpos)
  
  result = list(samp=samp, alg=alg, time=time, pos_summ = pos_summ)
  if (alg %in% c("ESSwS-MH", "ESSwS-ESS"))
  {
    result$betas = betas_out
  }
  
  return(result)
}

MH_betas = function(betas, lik_init, f, u)
{
  # Metropolis step for betas
  new_beta0 = rnorm(n = 1, mean = betas[1], sd = 0.1)
  new_beta1 = rnorm(n = 1, mean = betas[2], sd = 0.1)
  new_u = coal_samp_loglik(init = lik_init, f = f, beta0 = new_beta0, beta1 = new_beta1)
  if (new_u > u || log(runif(n = 1)) < new_u - u)
  {
    result = list(betas = c(new_beta0, new_beta1), u = new_u, ind = 1)
  }
  else
  {
    result = list(betas = betas, u = u, ind = 0)
  }
  return(result)
}

whiten_kappa = function(kappa, f, lik_init, cholC, invtcholC, loglikf, u, alpha, beta, prop_sd = 1)
{
  nu = (invtcholC * sqrt(kappa)) %*% f
  new_kappa = rlnorm(n = 1, meanlog = log(kappa), sdlog = prop_sd)
  new_f = (t(cholC) / sqrt(new_kappa)) %*% nu
  
  new_u = coal_loglik(init = lik_init, f = new_f) # plus other terms
  udiff = new_u - u
  priordiff = log_kappa_prior(kappa = new_kappa, alpha = alpha, beta = beta) -
    log_kappa_prior(kappa = kappa, alpha = alpha, beta = beta)
  propdiff = dlnorm(x = log(kappa), meanlog = log(new_kappa), sdlog = prop_sd) - 
    dlnorm(x = log(new_kappa), meanlog = log(kappa), sdlog = prop_sd)
  
  if (log(runif(n = 1)) < udiff + priordiff + propdiff)
  {
    result = list(kappa = new_kappa, f = new_f)
  }
  else
  {
    result = list(kappa = kappa, f = f)
  }
  
  return(result)
}

compute_pos_summ = function(samp_alg, loglikf, f, kappa, invC, alpha, beta,
                            betas=NULL, betas_prec=diag(1/100, 1/100))
{
  if (samp_alg == "ESS")
    loglik = loglikf(c(f, betas))
  else
    loglik = loglikf(f)
  
  
  logpri = log_mvnorm_prior(x = f, prec = invC * kappa) +
    log_kappa_prior(kappa = kappa, alpha = alpha, beta = beta)
  if (samp_alg %in% c("wS_MH", "wS_ESS"))
    logpri = logpri + log_betas_prior(betas, betas_prec)
  
  logpos = logpri + loglik
  
  return(list(loglik = loglik, logpri = logpri, logpos = logpos))
}

sampling_ESS = function(data, para, setting, init, samp_alg = "none", kappa_alg = "gibbs", surrogate = FALSE, verbose=TRUE)
{
  # pass the data and parameters
  lik_init = data$lik_init # f_offset = data$f_offset
  Ngrid = lik_init$ng+1
  alpha = para$alpha
  beta = para$beta
  invC = para$invC
  #rtEV = para$rtEV
  #EVC = para$EVC
  cholC = para$cholC
  
  # MCMC sampling setting
  #stepsz = setting$stepsz
  #Nleap  = setting$Nleap
  #if (alg=='aMALA')
  #  szkappa = setting$szkappa
  
  #if (alg=="HMC" | alg == "splitHMC")
  #{
  #  rand_leap = setting$rand_leap
  #}
  
  # storage of posterior samples
  NSAMP = setting$NSAMP
  #NBURNIN = setting$NBURNIN
  NBURNIN = 0
  
  #samp = matrix(NA, NSAMP - NBURNIN, Ngrid) # all parameters together
  #colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "kappa")
  
  # initialization
  #theta = init$theta
  f = init$theta[1:(Ngrid-1)]
  kappa = init$theta[Ngrid]
  u = init$u
  #du = init$du
  
  #acpi = 0
  acpt = rep(1, NSAMP - NBURNIN)
  
  fmat = matrix(NA, nrow = NSAMP - NBURNIN, ncol = length(f))
  kappas = rep(NA, NSAMP - NBURNIN)
  
  if (samp_alg %in% c("MH", "ESS"))
  {
    betas_out = matrix(NA, nrow = NSAMP - NBURNIN, ncol = 2)
  }
  
  # storage of log prior, log likelihood, and log posterior
  logpri = rep(NA, NSAMP - NBURNIN)
  loglik = rep(NA, NSAMP - NBURNIN)
  logpos = rep(NA, NSAMP - NBURNIN)
  
  if (samp_alg == "none")
  {
    ll = function(f) coal_loglik(init = lik_init, f = f)
  }
  else if (samp_alg == "fixed")
  {
    betas = para$betas
    ll = function(f) coal_samp_loglik(init = lik_init, f = f, beta0 = betas[1], beta1 = betas[2])
  }
  else if (samp_alg == "MH")
  {
    betas = init$betas
    ll = function(f) coal_samp_loglik(init = lik_init, f = f, beta0 = betas[1], beta1 = betas[2])
  }
  else if (samp_alg == "ESS")
  {
    betas = init$betas
    ll = function(f) coal_samp_loglik(init = lik_init, f = f[1:(Ngrid-1)], beta0 = f[Ngrid], beta1 = f[Ngrid+1])
  }
  
  if (kappa_alg == "whiten")
  {
    invtcholC = solve(t(cholC))
  }
  
  # start MCMC run
  start_time = Sys.time()
  cat('Running ESS', alg ,' sampling...\n')
  for(Iter in 1:NSAMP)
  {
    if (samp_alg == "none")
    {
      #ll = function(f) coal_loglik(init = lik_init, f = f)
      res = ESS(q_cur = f, l_cur = u, loglik = ll, cholC = cholC/sqrt(kappa))
      f = res$q
      
      #pos_summ = list(loglik = res$u)
      #pos_summ$logpri = log_mvnorm_prior(x = f, prec = invC * kappa) +
      #  log_kappa_prior(kappa = kappa, alpha = alpha, beta = beta)
      #pos_summ$logpos = pos_summ$logpri + pos_summ$loglik
      #res$pos_summ = pos_summ
    }
    else if (samp_alg == "fixed")
    {
      #ll = function(f) coal_samp_loglik(init = lik_init, f = f, beta0 = betas[1], beta1 = betas[2])
      res = ESS(q_cur = f, l_cur = u, loglik = ll, cholC = cholC/sqrt(kappa))
      f = res$q
      
      #pos_summ = list(loglik = res$u)
      #pos_summ$logpri = log_mvnorm_prior(x = f, prec = invC * kappa) +
      #  log_kappa_prior(kappa = kappa, alpha = alpha, beta = beta)
      #pos_summ$logpos = pos_summ$logpri + pos_summ$loglik
      #res$pos_summ = pos_summ
    }
    else if (samp_alg == "MH")
    {
      ## Metropolis step for betas
      #new_beta0 = rnorm(n = 1, mean = betas[1], sd = 0.1)
      #new_beta1 = rnorm(n = 1, mean = betas[2], sd = 0.1)
      #new_u = coal_samp_loglik(init = lik_init, f = f, beta0 = new_beta0, beta1 = new_beta1)
      #if (new_u > u || log(runif(n = 1)) < new_u - u)
      #{
      #  u = new_u
      #  betas = c(new_beta0, new_beta1)
      #}
      MH_betas(betas = betas, lik_init = lik_init, f = f, u = u)$betas
      
      #ll = function(f) coal_samp_loglik(init = lik_init, f = f, beta0 = betas[1], beta1 = betas[2])
      res = ESS(q_cur = f, l_cur = u, loglik = ll, cholC = cholC/sqrt(kappa))
      f = res$q
      
      #pos_summ = list(loglik = res$u)
      #pos_summ$logpri = log_mvnorm_prior(x = f, prec = invC * kappa) +
      #  log_kappa_prior(kappa = kappa, alpha = alpha, beta = beta)
      #pos_summ$logpos = pos_summ$logpri + pos_summ$loglik
      #res$pos_summ = pos_summ
    }
    else if (samp_alg == "ESS")
    {
      cholC_betas = matrix(0, nrow = nrow(cholC)+2, ncol = ncol(cholC)+2)
      cholC_betas[1:nrow(cholC), 1:ncol(cholC)] = cholC/sqrt(kappa)
      cholC_betas[nrow(cholC)+1, ncol(cholC)+1] = sqrt(100)
      cholC_betas[nrow(cholC)+2, ncol(cholC)+2] = sqrt(100)
      
      #ll = function(f) coal_samp_loglik(init = lik_init, f = f[1:(Ngrid-1)], beta0 = f[Ngrid], beta1 = f[Ngrid+1])
      res = ESS(q_cur = c(f, betas), l_cur = u, loglik = ll, cholC = cholC_betas)
      f = res$q[1:(Ngrid-1)]
      betas = res$q[Ngrid:(Ngrid+1)]
      
      #pos_summ = list(loglik = res$u)
      #pos_summ$logpri = log_mvnorm_prior(x = f, prec = invC * theta[Ngrid]) +
      #  log_mvnorm_prior(x = betas, prec = diag(c(1/100, 1/100))) +
      #  log_kappa_prior(kappa = kappa, alpha = alpha, beta = beta)
      #pos_summ$logpos = pos_summ$logpri + pos_summ$loglik
      #res$pos_summ = pos_summ
    }
    else
    {
      stop('The ESS subalgorithm is not in the list!')
    }
    
    #acpi <- acpi+res$Ind
    
    if (kappa_alg == "gibbs")
    {
      kappa <- rgamma(1,alpha+(Ngrid-1)/2, beta + t(f) %*% invC %*% f/2)
    }
    else if (kappa_alg == "whiten")
    {
      kappa_res <- whiten_kappa(kappa = kappa, f = f, lik_init = lik_init,
                                cholC = cholC, invtcholC = invtcholC,
                                loglikf = ll, u = u, alpha = alpha, beta = beta)
      kappa <- kappa_res$kappa
      f <- kappa_res$f
    }
    else
      stop("Kappa operator not recognized.")
    
    pos_summ = compute_pos_summ(samp_alg = samp_alg, loglikf = ll, f = f, kappa = kappa,
                                invC = invC, alpha = alpha, beta = beta, betas = betas)
    u <- pos_summ$loglik
    
    # save posterior samples after burnin
    if (Iter > NBURNIN)
    {
      fmat[Iter - NBURNIN, ] <- f
      kappas[Iter - NBURNIN] <- kappa
      
      if (samp_alg %in% c("MH", "ESS"))
      {
        betas_out[Iter - NBURNIN, ] <- betas
      }
      #acpt[Iter - NBURNIN] <- res$Ind
      
      logpri[Iter - NBURNIN] <- pos_summ$logpri
      loglik[Iter - NBURNIN] <- pos_summ$loglik
      logpos[Iter - NBURNIN] <- pos_summ$logpos
      
      #if (samp_alg %in% c("MH", "ESS"))
      #{
      #  betas_out = rbind(betas_out, betas, deparse.level = 0)
      #}
    }
    
    if (verbose && Iter %% 100 == 0)
    {
      cat(Iter, ' iterations have been finished!\n' )
      cat('Online acceptance rate is ', 1,'\n')
      #acpi=0
    }
  }
  stop_time <- Sys.time()
  time <- stop_time-start_time
  cat('\nTime consumed : ',time)
  #acpt <- acpt/(NSAMP-NBURNIN)
  cat('\nFinal Acceptance Rate: ', sum(acpt) / (NSAMP-NBURNIN),'\n')
  
  pos_summ = data.frame(acpt=acpt, logpri = logpri, loglik = loglik, logpos = logpos)
  
  if (samp_alg %in% c("MH", "ESS"))
  {
    samp = cbind(fmat, kappas, betas_out)
    colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "kappa", "beta0", "beta1")
  }
  else
  {
    samp = cbind(fmat, kappas)
    colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "kappa")
  }
  
  return(list(samp=samp, alg="ESS", time=time, pos_summ = pos_summ, samp_alg = samp_alg))
}

burnin_subsample = function(res, burnin = 0, subsample = 1)
{
  nsamp = dim(res$samp)[1]
  indices = seq(burnin+1, nsamp, by=subsample)
  samp = res$samp[indices, ]
  pos_summ = res$pos_summ[indices, ]
  
  cleaned_result = list(samp=samp, pos_summ = pos_summ, grid = res$grid)
}

calculate_estimates = function(logfmat, params, grid)
{
  logfmed = apply(logfmat, MARGIN = 2, median)
  logflow = apply(logfmat, MARGIN = 2, function(x) quantile(x, .025))
  logfhi  = apply(logfmat, MARGIN = 2, function(x) quantile(x, .975))
  
  logfmed_fun = stepfun(grid, c(0, logfmed, 0))
  logflow_fun = stepfun(grid, c(0, logflow, 0))
  logfhi_fun  = stepfun(grid, c(0, logfhi,  0))
  
  fmed = apply(exp(logfmat), MARGIN = 2, median)
  flow = apply(exp(logfmat), MARGIN = 2, function(x) quantile(x, .025))
  fhi  = apply(exp(logfmat), MARGIN = 2, function(x) quantile(x, .975))
  
  fmed_fun = stepfun(grid, c(0, fmed, 0))
  flow_fun = stepfun(grid, c(0, flow, 0))
  fhi_fun  = stepfun(grid, c(0, fhi,  0))
  
  pmed = apply(params, MARGIN = 2, median)
  plow = apply(params, MARGIN = 2, function(x) quantile(x, .025))
  phi  = apply(params, MARGIN = 2, function(x) quantile(x, .975))
  
  return(list(logfmed = logfmed, logflow = logflow, logfhi = logfhi,
              logfmed_fun = logfmed_fun, logflow_fun = logflow_fun,
              logfhi_fun = logfhi_fun, fmed = fmed, flow = flow, fhi = fhi,
              fmed_fun = fmed_fun, flow_fun = flow_fun, fhi_fun = fhi_fun,
              pmed = pmed, plow = plow, phi = phi))
}

# wrapper that encapsulates sampler above with good defaults
#' @export
mcmc_sampling = function(data, alg, nsamp, nburnin, Ngrid=100, nugget="1,1", prec_alpha = 1e-2, prec_beta = 1e-2,
                         TrjL=NULL, Nleap=NULL, szkappa=NULL, rand_leap=NULL, betas=c(0,0), samp_alg = "none", kappa_alg = "gibbs")
{
  # add ability to parse genealogy objects as well as lists
  samp_times = data$samp_times
  n_sampled  = data$n_sampled
  coal_times = data$coal_times
  
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
  
  grid = seq(grid_bds[1],grid_bds[2],length.out=Ngrid)
  intl = grid[2]-grid[1]
  midpts = grid[-1]-intl/2
  
  # initialize likelihood calculation
  lik_init = coal_lik_init(samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times, grid=grid)
  
  # calculate intrinsic precision matrix
  invC <- Q_matrix(as.matrix(midpts),0,1)
  
  # fudge to be able to compute the cholC
  if (nugget == "1,1")
    invC[1,1] <- invC[1,1]+.0001 # nugget at (1,1)
  else if (nugget == "diag")
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
  theta = rep(1,Ngrid)
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
    u  = coal_loglik(init = lik_init, f = theta[-Ngrid])
    du = NULL
  }
  else if (alg == "ESSwS")
  {
    u  = coal_samp_loglik(init = lik_init, f = theta[-Ngrid], beta0 = betas[1], beta1 = betas[2])
    du = NULL
  }
  else if (alg == "ESSwS-MH")
  {
    u  = coal_samp_loglik(init = lik_init, f = theta[-Ngrid], beta0 = betas[1], beta1 = betas[2])
    du = NULL
  }
  else if (alg == "ESSwS-ESS")
  {
    u  = coal_samp_loglik(init = lik_init, f = theta[-Ngrid], beta0 = betas[1], beta1 = betas[2])
    du = NULL
  }
  else
  {
    stop('The algorithm is not in the list!')
  }
  
  # MCMC sampling preparation
  data = list(lik_init=lik_init)
  para = list(alpha=prec_alpha,beta=prec_beta,invC=invC,rtEV=rtEV,EVC=EVC,cholC=cholC, betas=betas)
  setting = list(stepsz=stepsz,Nleap=Nleap,NSAMP=nsamp,NBURNIN=nburnin,szkappa = szkappa,rand_leap=rand_leap)
  init = list(theta=theta,u=u,du=du, betas=betas)
  
  # Run MCMC sampler
  if (alg == "ESS")
  {
    res_MCMC = sampling_ESS(data = data, para = para, setting = setting,
                            init = init, samp_alg = samp_alg, kappa_alg = kappa_alg)
  }
  else
  {
    res_MCMC = sampling(data = data, para = para, alg = alg, setting = setting,
                        init = init)
  }
  
  cleaned_res = burnin_subsample(res = res_MCMC, burnin = nburnin)
  
  logfmat = cleaned_res$samp[,1:(Ngrid-1)]
  if (alg == "ESS" && samp_alg %in% c("MH", "ESS"))
  {
    params = cleaned_res$samp[,Ngrid:(Ngrid+2)]
  }
  else
  {
    params = matrix(cleaned_res$samp[,Ngrid])
  }
  estimates = calculate_estimates(logfmat = logfmat, params = params, grid = grid)
  
  res_MCMC$cleaned_res = cleaned_res
  res_MCMC$estimates = estimates
  
  res_MCMC$med = estimates$fmed
  res_MCMC$low = estimates$flow
  res_MCMC$hi = estimates$fhi
  
  res_MCMC$med_fun = estimates$fmed_fun
  res_MCMC$low_fun = estimates$flow_fun
  res_MCMC$hi_fun = estimates$fhi_fun
  
  return(res_MCMC)
}

sampling_old = function(data, para, alg, setting, init, verbose=TRUE)
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
  
  
  # storage of posterior samples
  NSAMP = setting$NSAMP
  NBURNIN = setting$NBURNIN
  samp = matrix(NA,NSAMP-NBURNIN,Ngrid) # all parameters together
  acpi = 0
  acpt = 0
  
  # initialization
  theta = init$theta
  u = init$u
  du = init$du
  
  # start MCMC run
  start_time = Sys.time()
  cat('Running ', alg ,' sampling...\n')
  for(Iter in 1:NSAMP)
  {
    if(verbose&&Iter%%100==0)
    {
      cat(Iter, ' iterations have been finished!\n' )
      cat('Online acceptance rate is ',acpi/100,'\n')
      acpi=0
    }
    
    # sample the whole parameter
    #tryCatch({res=switch(alg,
    #                     HMC=eval(parse(text='HMC'))(theta,u,du,function(theta,grad=FALSE)U(theta,lik_init,invC,alpha,beta,grad),stepsz,Nleap,rand_leap),
    #                     splitHMC=eval(parse(text='splitHMC'))(theta,u,du,function(theta,grad=FALSE)U_split(theta,lik_init,invC,alpha,beta,grad),rtEV,EVC,stepsz,Nleap,rand_leap),
    #                     MALA=eval(parse(text='MALA'))(theta,u,du,function(theta,grad=FALSE)U(theta,lik_init,invC,alpha,beta,grad),stepsz),
    #                     aMALA=eval(parse(text='aMALA'))(theta,u,function(theta,grad=FALSE)U_kappa(theta,lik_init,invC,alpha,beta,grad),function(theta)Met(theta,lik_init,invC),szkappa,stepsz),
    #                     ESS=eval(parse(text='ESS'))(theta[-Ngrid],u,function(f)coal_loglik(lik_init,f),cholC/sqrt(theta[Ngrid])),
    #                     stop('The algorithm is not in the list!'));
    #          theta[1:(Ngrid-(alg=='ESS'))]=res$q;u=res[[2]];if(any(grepl(alg,c('HMC','splitHMC','MALA'))))du=res$du;
    #          acpi=acpi+res$Ind}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    res=switch(alg,
               HMC=eval(parse(text='HMC'))(theta,u,du,function(theta,grad=FALSE)U(theta,lik_init,invC,alpha,beta,grad),stepsz,Nleap,rand_leap),
               splitHMC=eval(parse(text='splitHMC'))(theta,u,du,function(theta,grad=FALSE)U_split(theta,lik_init,invC,alpha,beta,grad),rtEV,EVC,stepsz,Nleap,rand_leap),
               MALA=eval(parse(text='MALA'))(theta,u,du,function(theta,grad=FALSE)U(theta,lik_init,invC,alpha,beta,grad),stepsz),
               aMALA=eval(parse(text='aMALA'))(theta,u,function(theta,grad=FALSE)U_kappa(theta,lik_init,invC,alpha,beta,grad),function(theta)Met(theta,lik_init,invC),szkappa,stepsz),
               ESS=eval(parse(text='ESS'))(theta[-Ngrid],u,function(f)coal_loglik(lik_init,f),cholC/sqrt(theta[Ngrid])),
               stop('The algorithm is not in the list!'));
    theta[1:(Ngrid-(alg=='ESS'))]=res$q
    u=res[[2]]
    if(any(grepl(alg,c('HMC','splitHMC','MALA'))))du=res$du
    acpi=acpi+res$Ind
    # Gibbs sample kappa for ESS
    if(alg=='ESS')
      theta[Ngrid]=rgamma(1,alpha+(Ngrid-1)/2,beta+t(theta[-Ngrid])%*%invC%*%theta[-Ngrid]/2)
    
    # save posterior samples after burnin
    if(Iter>NBURNIN)
    {
      samp[Iter-NBURNIN,]<-theta
      acpt<-acpt+res$Ind
    }
    
  }
  stop_time = Sys.time()
  time = stop_time-start_time
  cat('\nTime consumed : ',time)
  acpt = acpt/(NSAMP-NBURNIN)
  cat('\nFinal Acceptance Rate: ',acpt,'\n')
  
  return(list(samp=samp,time=time,acpt=acpt))
}
