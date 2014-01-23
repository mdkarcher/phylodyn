coal_lik_init = function(samp_times, n_sampled, coal_times, grid)
{
  ns = length(samp_times)
  nc = length(coal_times)
  ng = length(grid)-1
  
  if (length(samp_times) != length(n_sampled))
    stop("samp_times vector of differing length than n_sampled vector.")
  
  if (length(coal_times) != sum(n_sampled) - 1)
    stop("Incorrect length of coal_times: should be n_sampled - 1.")
  
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
    gridrep[i] = sum(t >= grid[i] & t < grid[i+1])
  
  C = 0.5 * l * (l-1)
  D = diff(t)
  
  y = rep(0, length(D))
  y[t[-1] %in% coal_times] = 1
  
  return(list(t=t, l=l, C=C, D=D, y=y, gridrep=gridrep, ng=ng, args=list(samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times, grid=grid)))
}

coal_lik_eval = function(init, N, log=TRUE)
{
  if (init$ng != length(N))
    stop(paste("Incorrect length for N; should be", init$ng))
  
  N = rep(N, init$gridrep)
  
  rat = init$C / N
  #print(rat)
  
  lls = init$y * log(rat) - init$D * rat
  #print(lls)
  
  ll = sum(lls[!is.nan(lls)])
  
  if (log)
    result = sum(ll)
  else
    result = exp(sum(ll))
  
  return(result)
}