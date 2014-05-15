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

traj_beta = function(t, traj, beta)
{
  return(traj(t)^beta)
}

unif_traj = function(t, level=100)
{
  n = length(t)
  return(rep(level,n))
}

unif_traj_inv = function(t, level=100)
{
  n = length(t)
  return(rep(1/level,n))
}

exp_traj = function(t, scale=1000)
{
  return(scale * exp(-t))
}

exp_traj_inv = function(t, scale=1000)
{
  return(1/(scale * exp(-t)))
}

boombust_traj = function(t, bust=1, scale=1000)
{
  #bust = 1
  result = rep(0, length(t))
  result[t <= bust] = scale*exp(t[t <= bust]-bust)
  result[t >  bust] = scale*exp(bust-t[t >  bust])
  return(result)
}

boombust_traj_inv = function(t, bust=1, scale=1000)
{
  return(1/boombust_traj(t, bust, scale))
}

cyclic_traj = function(t)
{
  result = rep(0, length(t))
  result[(t %% 10) <= 5] = 200*exp(-(t[(t %% 10) <= 5] %% 10) / 2)
  result[(t %% 10) >  5] = 200*exp((t[(t %% 10) >  5] %% 10) / 2 - 5)
  return(result)
}

cyclic_traj_inv = function(t)
{
  return(1/cyclic_traj(t))
}

steep_cyc_traj = function(t)
{
  result = rep(0, length(t))
  result[(t %% 10) <= 5] = 20 + 1980*exp(-(t[(t %% 10) <= 5] %% 10) * 5)
  result[(t %% 10) >  5] = 20 + 1980*exp(((t[(t %% 10) >  5] %% 10) - 10) * 5)
  return(result)
}

steep_cyc_traj_inv = function(t)
{
  return(1/steep_cyc_traj(t))
}

sloped_traj = function(t)
{
  result = rep(0, length(t))
  result[(t %% 10) <= 5] = 100 + 900*exp(-(t[(t %% 10) <= 5] %% 10) * 1)
  result[(t %% 10) >  5] = 100 + 900*exp(((t[(t %% 10) >  5] %% 10) - 10) * 1)
  return(result)
}

sloped_traj_inv = function(t)
{
  return(1/sloped_traj(t))
}

mesa_traj = function(t, a=2, b=3)
{
  result = rep(0, length(t))
  result[(t %% 10) <= a | (t %% 10) > 10-a] = 1000
  result[(t %% 10) >  a & (t %% 10) <= 5] = 100 + 900*exp((a-((t[(t %% 10) >  a & (t %% 10) <= 5] %% 10))) * b)
  result[(t %% 10) >  5 & (t %% 10) <= 10-a] = 100 + 900*exp(((t[(t %% 10) >  5 & (t %% 10) <= 10-a] %% 10) - 10 + a) * b)
  return(result)
}

mesa_traj_inv = function(t, a=2, b=3)
{
  return(1/mesa_traj(t, a, b))
}

logistic_traj = function(t, offset=0, a=2)
{
  t = t + offset
  result = rep(0, length(t))
  result[(t %% 12) <= 6] = 10 + 90/(1+exp((3-(t[(t %% 12) <= 6] %% 12)) * a))
  result[(t %% 12) >  6] = 10 + 90/(1+exp(((t[(t %% 12) >  6] %% 12) - 12 + 3) * a))
  return(result)
}

logistic_traj_inv = function(t, offset=0, a=2)
{
  return(1/logistic_traj(t, offset, a))
}

# sre = function(hat_f, true_f, s)
# {
#   return(sum( abs(hat_f(s)-true_f(s)) / true_f(s) ))
# }
# 
# mre = function(hat_f, true_f, s)
# {
#   return(sum( (hat_f(s)-true_f(s)) / true_f(s) ))
# }
# 
# mrw = function(hat, hat975, hat025, s)
# {
#   K = length(s)
#   return(sum( abs(hat975(s) - hat025(s)) / (K * hat(s)) ))
# }
# 
# envelope = function(true_f, hat975, hat025, s)
# {
#   K = length(s)
#   return(sum(true_f(s) < hat975(s) & true_f(s) > hat025(s)) / K)
# }
# 
# variation = function(hat_f, s)
# {
#   y = hat_f(s)
#   return(sum(abs(diff(y))))
# }