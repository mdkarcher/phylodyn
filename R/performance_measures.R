bias = function(INLA_out, traj, hilim=Inf, lolim=0, yhilim=Inf, ylolim=0)
{
  mod = INLA_out$result$summary.random$time
  mask = mod$ID <= hilim & mod$ID >= lolim & traj(mod$ID) <= yhilim & traj(mod$ID) >= lolim
  grid_pts = mod$ID[mask]
  n = length(grid_pts)
  med = exp(-mod$"0.5quant"[mask])
  truth = traj(grid_pts)
  result = sum( (med - truth)/truth )
  return(list(tot = result, avg = result / n ))
}

dev = function(INLA_out, traj, hilim=Inf, lolim=0, yhilim=Inf, ylolim=0)
{
  mod = INLA_out$result$summary.random$time
  mask = mod$ID <= hilim & mod$ID >= lolim & traj(mod$ID) <= yhilim & traj(mod$ID) >= lolim
  grid_pts = mod$ID[mask]
  n = length(grid_pts)
  med = exp(-mod$"0.5quant"[mask])
  truth = traj(grid_pts)
  result = sum( abs(med - truth)/truth )
  return(list(tot = result, avg = result / n ))
}

relwid = function(INLA_out, traj, hilim=Inf, lolim=0, yhilim=Inf, ylolim=0)
{
  mod = INLA_out$result$summary.random$time
  mask = mod$ID <= hilim & mod$ID >= lolim & traj(mod$ID) <= yhilim & traj(mod$ID) >= lolim
  grid_pts = mod$ID[mask]
  n = length(grid_pts)
  lo = exp(-mod$"0.975quant"[mask])
  hi = exp(-mod$"0.025quant"[mask])
  truth = traj(grid_pts)
  result = sum( (hi - lo)/truth )
  return(list(tot = result, avg = result / n ))
}

relwid_emp = function(INLA_out, hilim=Inf, lolim=0)
{
  mod = INLA_out$result$summary.random$time
  mask = mod$ID <= hilim & mod$ID >= lolim
  grid_pts = mod$ID[mask]
  n = length(grid_pts)
  med = exp(-mod$"0.5quant"[mask])
  lo = exp(-mod$"0.975quant"[mask])
  hi = exp(-mod$"0.025quant"[mask])
  result = sum( (hi - lo)/med )
  return(list(tot = result, avg = result / n ))
}

envelope = function(INLA_out, traj, hilim=Inf, lolim=0, yhilim=Inf, ylolim=0)
{
  mod = INLA_out$result$summary.random$time
  mask = mod$ID <= hilim & mod$ID >= lolim & traj(mod$ID) <= yhilim & traj(mod$ID) >= lolim
  grid_pts = mod$ID[mask]
  n = length(grid_pts)
  lo = exp(-mod$"0.975quant"[mask])
  hi = exp(-mod$"0.025quant"[mask])
  truth = traj(grid_pts)
  result = sum(truth < hi & truth > lo)
  return(list(tot = result, avg = result / n))
}

variation = function(INLA_out, hilim=Inf, lolim=0)
{
  mod = INLA_out$result$summary.random$time
  mask = mod$ID <= hilim & mod$ID >= lolim
  grid_pts = mod$ID[mask]
  n = length(grid_pts)
  med = exp(-mod$"0.5quant"[mask])
  result = sum(abs(diff(med)) / med[1:(n-1)])
  return(list(tot=result, avg=result/n))
}
