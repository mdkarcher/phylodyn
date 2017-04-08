#### Coalescent likelihood functions ####
# Compute log likelihood of coalescent model, energy function in HMC algorithms, and metric tensor needed in aMALA.

coal_lik_init = function(samp_times, n_sampled, coal_times, grid)
{
  ns = length(samp_times)
  nc = length(coal_times)
  ng = length(grid)-1
  
  if (length(samp_times) != length(n_sampled))
    stop("samp_times vector of differing length than n_sampled vector.")
  
  if (length(coal_times) != sum(n_sampled) - 1)
    stop("Incorrect length of coal_times: should be sum(n_sampled) - 1.")
  
  if (max(samp_times, coal_times) > max(grid))
    stop("Grid does not envelop all sampling and/or coalescent times.")
  
  t = sort(unique(c(samp_times, coal_times, grid)))
  l = rep(0, length(t))
  
  for (i in 1:ns)
    l[t >= samp_times[i]] = l[t >= samp_times[i]] + n_sampled[i]
  
  for (i in 1:nc)
    l[t >= coal_times[i]] = l[t >= coal_times[i]] - 1
  
  #print(l)
  
  if (sum((l < 1) & (t >= min(samp_times))) > 0)
    stop("Number of active lineages falls below 1 after the first sampling point.")
  
  mask = l > 0
  t = t[mask]
  l = utils::head(l[mask], -1)
  
  gridrep = rep(0, ng)
  for (i in 1:ng)
    gridrep[i] = sum(t > grid[i] & t <= grid[i+1])
  
  C = 0.5 * l * (l-1)
  D = diff(t)
  
  y = rep(0, length(D))
  y[t[-1] %in% coal_times] = 1
  
  bins = cut(x = samp_times, breaks = t,
                include.lowest = TRUE)
  tab <- stats::aggregate(n_sampled ~ bins, FUN = sum, labels = FALSE)
  count <- rep(0, length(D))
  count[as.numeric(tab$bins)] <- tab$n_sampled
  count[utils::head(t, -1) >= max(samp_times)] <- NA
  
  rep_idx = cumsum(gridrep)
  rep_idx = cbind(rep_idx-gridrep+1,rep_idx)
  
  return(list(t=t, l=l, C=C, D=D, y=y, count=count, gridrep=gridrep, ns=sum(n_sampled), nc=nc, ng=ng, rep_idx=rep_idx, args=list(samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times, grid=grid)))
}

coal_loglik = function(init, f, grad=FALSE)
{
  if (init$ng != length(f))
    stop(paste("Incorrect length for f; should be", init$ng))
  
  f = rep(f, init$gridrep)
  
  llnocoal = init$D * init$C * exp(-f)
  
  if (!grad)
  {  
    lls = -init$y * f - llnocoal
    #print(lls)
    
    ll = sum(lls[!is.nan(lls)])
    
    return(ll)
  }
  else
  {  
    dll = apply(init$rep_idx,1,function(idx)sum(-init$y[idx[1]:idx[2]]+llnocoal[idx[1]:idx[2]])) # gradient of log-likelihood wrt f_midpts
    
    return(dll)
  }
}

samp_loglik = function(init, fs, betas)
{
  fs = as.matrix(fs)
  
  if (init$ng != dim(fs)[1])
    stop(paste("Incorrect number of rows for fs; should be", init$ng))
  
  if (length(betas) != dim(fs)[2] + 1)
    stop(paste("Incompatible number of betas and columns for fs"))
  
  beta0 = betas[1]
  betas = utils::tail(betas, -1)
  
  #f = rep(f, init$gridrep)
  fs = as.matrix(fs[rep(1:init$ng, init$gridrep),])
  
  #llsampevents = beta1 * init$count * f
  #llsampnoevents = init$D * exp(beta0) * exp(f)^beta1
  #llsamp = init$ns * beta0 + sum(llsampevents[!is.na(init$count)]) - sum(llsampnoevents[!is.na(init$count)])
  
  fs_betas = fs %*% betas
  
  llsampevents = init$count * fs_betas
  llsampnoevents = init$D * exp(beta0 + fs_betas)
  llsamp = init$ns * beta0 + sum(llsampevents[!is.na(init$count)]) - sum(llsampnoevents[!is.na(init$count)])
  #print(init$ns * beta0)
  #print(sum(llsampevents[!is.na(init$count)]))
  #print(-1 * sum(llsampnoevents[!is.na(init$count)]))
  
  return(llsamp)
}

coal_samp_loglik = function(init, f, beta0, beta1)
{
  if (init$ng != length(f))
    stop(paste("Incorrect length for f; should be", init$ng))
  
  f = rep(f, init$gridrep)
  
  llnocoal = init$D * init$C * exp(-f)
  
  lls = -init$y * f - llnocoal
  
  llsampevents = beta1 * init$count * f
  llsampnoevents = init$D * exp(beta0) * exp(f)^beta1
  llsamp = init$ns * beta0 + sum(llsampevents[!is.na(init$count)]) - sum(llsampnoevents[!is.na(init$count)])
  
  llcoal = sum(lls[!is.nan(lls)])
  
  return(llcoal + llsamp)
}

coal_samp_fns_loglik = function(init, f, fs, beta0, beta1, betas)
{
  if (init$ng != length(f))
    stop(paste("Incorrect length for fs; should be", init$ng))
  
  if (init$ng != dim(fs)[1])
    stop(paste("Incorrect number of rows for fs; should be", init$ng))
  
  f = rep(f, init$gridrep)
  fs = fs[rep(1:init$ng, init$gridrep),]
  
  llnocoal = init$D * init$C * exp(-f)
  
  lls = -init$y * f - llnocoal
  
  llsampevents = beta1 * init$count * f + diag(init$count) %*% fs %*% betas
  llsampnoevents = init$D * exp(beta0 + f * beta1 + fs %*% betas)
  llsamp = init$ns * beta0 + sum(llsampevents[!is.na(init$count)]) - sum(llsampnoevents[!is.na(init$count)])
  
  llcoal = sum(lls[!is.nan(lls)])
  
  return(llcoal + llsamp)
}

U = function(theta, init, invC, alpha, beta, grad=FALSE)
{
  D = length(theta)
  f = theta[-D]
  tau = theta[D]
  invCf = invC %*% f
  if(!grad)
  {
    loglik = coal_loglik(init, f)
    logpri = ((D-1)/2+alpha)*tau - (t(f)%*%invCf/2+beta)*exp(tau)
    return(list(loglik = -loglik, logpri = -logpri, logpos = -(loglik+logpri)))
  }
  else
  {
    dloglik = c(coal_loglik(init, f, grad),0)
    dlogpri = c(-invCf*exp(tau),((D-1)/2+alpha)-(t(f)%*%invCf/2+beta)*exp(tau))
    return(list(dloglik = -dloglik, dlogpri = -dlogpri, dlogpos = -(dloglik+dlogpri)))
  }
}

log_f_prior_kappa = function(f, kappa, invC, alpha, beta)
{
  D = length(f) + 1
  invCf = invC %*% f
  
  return(((D-1)/2+alpha-1)*log(kappa) - (t(f)%*%invCf/2+beta)*kappa)
}

U_kappa = function(theta, init, invC, alpha, beta, grad=FALSE)
{
  D = length(theta)
  f = theta[-D]
  kappa=theta[D]
  invCf = invC %*% f
  if(!grad)
  {
    loglik = coal_loglik(init, f)
    logpri = ((D-1)/2+alpha-1)*log(kappa) - (t(f)%*%invCf/2+beta)*kappa
    return(list(loglik = -loglik, logpri = -logpri, logpos = -(loglik+logpri)))
  }
  else
  {
    dloglik = c(coal_loglik(init, f, grad),0)
    dlogpri = c(-invCf*kappa,((D-1)/2+alpha-1)/kappa-(t(f)%*%invCf/2+beta))
    return(list(dloglik = -dloglik, dlogpri = -dlogpri, dlogpos = -(dloglik+dlogpri)))
  }
}

U_split = function(theta, init, invC, alpha, beta, grad=FALSE)
{
  D=length(theta)
  f=theta[-D]
  tau=theta[D]
  invCf=invC%*%f
  if(!grad)
  {
    loglik = coal_loglik(init, f)
    logpri = ((D-1)/2+alpha)*tau - (t(f)%*%invCf/2+beta)*exp(tau)
    return(list(loglik = -loglik, logpri = -logpri, logpos = -(loglik+logpri)))
  }
  else
  {
    dU_res = -c(coal_loglik(init, f, grad),((D-1)/2+alpha)-beta*exp(tau))
    return(dU_res)	
  }
}

precBM = function(times, delta=1e-6)
{
  D = length(times)
  
  diff1 <- diff(times)
  diff1[diff1==0] <- delta
  
  diff <- 1/diff1
  Q <- spam::spam(0,D,D)
  if (D>2)
    Q[cbind(1:D,1:D)] <- c(diff[1]+ifelse(times[1]==0,1/delta,1/times[1]),diff[1:(D-2)]+diff[2:(D-1)],diff[D-1])
  else
    Q[cbind(1:D,1:D)] <- c(diff[1]+ifelse(times[1]==0,1/delta,1/times[1]),diff[D-1])
  
  Q[cbind(1:(D-1),2:D)] = -diff[1:(D-1)]; Q[cbind(2:D,1:(D-1))]=-diff[1:(D-1)]
  return(Q)
}

Met = function(theta, init, invC)
{
  D = length(theta)
  f = theta[-D]
  kappa=theta[D]
  
  f = rep(f, init$gridrep)
  llnocoal = init$D * init$C * exp(-f)
  diagF = apply(init$rep_idx,1,function(idx)sum(llnocoal[idx[1]:idx[2]]))
  
  G = invC*kappa + diag(diagF)
  return(G)
}
