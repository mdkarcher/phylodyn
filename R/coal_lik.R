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
  l = head(l[mask], -1)
  
  gridrep = rep(0, ng)
  for (i in 1:ng)
    gridrep[i] = sum(t > grid[i] & t <= grid[i+1])
  
  C = 0.5 * l * (l-1)
  D = diff(t)
  
  y = rep(0, length(D))
  y[t[-1] %in% coal_times] = 1
  
  rep_idx = cumsum(gridrep)
  rep_idx = cbind(rep_idx-gridrep+1,rep_idx)
  
  return(list(t=t, l=l, C=C, D=D, y=y, gridrep=gridrep, ng=ng, rep_idx=rep_idx, args=list(samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times, grid=grid)))
}

coal_loglik = function(init, f, grad=FALSE)
{
  if (init$ng != length(f))
    stop(paste("Incorrect length for f; should be", init$ng))
  
  f = rep(f, init$gridrep)
  
  llnocoal = init$D * init$C * exp(-f)
  
  if (!grad)
  {  
    lls = - init$y * f - llnocoal
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
    return(-(loglik+logpri))
  }
  else
  {
    dloglik = c(coal_loglik(init, f, grad),0)
    dlogpri = c(-invCf*exp(tau),((D-1)/2+alpha)-(t(f)%*%invCf/2+beta)*exp(tau))
    return(-(dloglik+dlogpri))
  }
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
    return(-(loglik+logpri))
  }
  else
  {
    dloglik = c(coal_loglik(init, f, grad),0)
    dlogpri = c(-invCf*kappa,((D-1)/2+alpha-1)/kappa-(t(f)%*%invCf/2+beta))
    return(-(dloglik+dlogpri))
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
    return(-(loglik+logpri))
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
  Q<-spam(0,D,D)
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
