pref_sample = function(f, c=1, lim=c(0,1), grid.len = 1000, beta=1)
{
  grid = seq(lim[1], lim[2], length.out=grid.len)
  eval = f(grid)^beta
  upper = max(eval)
  
  nfull = rpois(1, c * upper * diff(lim))
  pfull = sort(runif(nfull, lim[1], lim[2]))
  
  pthin = NULL
  for (p in pfull)
  {
    if (runif(1) < f(p)^beta / upper)
      pthin = c(pthin, p)
  }
  #nthin = length(pthin)
  
  return(pthin)
}
