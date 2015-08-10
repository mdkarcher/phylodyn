shade_band = function(x, ylo, yhi, xlim=NULL, col="gray")
{
  if (is.null(xlim))
    xlim = c(0, Inf)
  mask = x >= min(xlim) & x <= max(xlim)
  
  x = x[mask]
  ylo = ylo[mask]
  yhi = yhi[mask]
  
  polygon(c(x, rev(x)), c(yhi, rev(ylo)), col=col, border=NA)
}

plot_BNPR_exper = function(INLA_out, traj=NULL, args=NULL, xlim=NULL, ylim=NULL,
                           nbreaks=40, lty=1, lwd=2, col="black",
                           xlab="Time", ylab="Effective Population Size",
                           subdir="time", label_events=FALSE, ...)
{
  mod = INLA_out$result$summary.random[[subdir]]
  
  grid = mod$"ID"
  if (is.null(xlim))
  {
    xlim=c(max(grid),0)
  }
  mask = grid >= min(xlim) & grid <= max(xlim)
  
  if (is.null(ylim))
  {
    ymax=max(exp(-mod$"0.025quant"[mask]))
    ymin=min(exp(-mod$"0.975quant"[mask]))
  }
  else
  {
    ymin=min(ylim)
    ymax=max(ylim)
  }
  
  if (is.null(args))
    ylim=c(ymin, ymax)
  else
  {
    yspan=ymax/ymin
    yextra=yspan^(1/10)
    ylim=c(ymin/(yextra^1.25), ymax)
  }
  plot(1,1,type="n",log="y",
       xlab=xlab, ylab=ylab,
       xlim=xlim, ylim=ylim, ...)
  #lines(grid,exp(-mod$"0.975quant"),lwd=2.5,col="blue",lty=2)
  #lines(grid,exp(-mod$"0.025quant"),lwd=2.5,col="blue",lty=2)
  shade_band(grid[mask], exp(-mod$"0.975quant")[mask],exp(-mod$"0.025quant")[mask], col="lightgray")
  lines(grid[mask], exp(-mod$"0.5quant")[mask], lwd=lwd, col=col, lty=lty)
  if (!is.null(traj))
    lines(grid[mask], traj(grid)[mask], lwd=2, lty=3)
  if (!is.null(args))
  {
    samps=args$s[args$event==0]
    samps = samps[samps < max(xlim) & samps > min(xlim)]
    coals=args$s[args$event==1]
    coals = coals[coals < max(xlim) & coals > min(xlim)]
    
    breaks = seq(min(xlim), max(xlim), length.out=nbreaks)
    h_samp = hist(samps, breaks=breaks, plot=FALSE)
    h_coal = hist(coals, breaks=breaks, plot=FALSE)
    
    hist2heat(h_samp, y=ymin/yextra^0.5, wd=10)
    hist2heat(h_coal, y=ymin/yextra, wd=10)
    
    #points(samps, rep(ymin/yextra^0.5, length(samps)), pch=3) 
    #points(coals, rep(ymin/yextra, length(coals)), pch=4)
    
    if (label_events)
    {
      text(x = xlim[2], y = ymin/(yextra^0.20), labels = "Sampling events", adj = c(1, 0), cex = 0.7)
      text(x = xlim[2], y = ymin/(yextra^1.25), labels = "Coalescent events", adj = c(1, 1), cex = 0.7)
    }
  }
}

plot_BNPR = function(BNPR_out, traj=NULL, xlim=NULL, ...)
{
  mod = BNPR_out$result$summary.random$time
  
  grid = mod$"ID"
  if (is.null(xlim))
  {
    xlim=c(max(grid),0)
  }
  plot(grid,exp(-mod$"0.5quant"),type="l",lwd=2.5,col="blue",log="y",
       xlab="Time (past to present)",ylab="Scaled Effective Pop. Size",
       xlim=xlim,ylim=c(min(exp(-mod$"0.975quant"[grid > min(xlim) & grid < max(xlim)])),
                        max(exp(-mod$"0.025quant"[grid > min(xlim) & grid < max(xlim)]))), ...)
  lines(grid,exp(-mod$"0.975quant"),lwd=2.5,col="blue",lty=2)
  lines(grid,exp(-mod$"0.025quant"),lwd=2.5,col="blue",lty=2)
  if (!is.null(traj))
    lines(grid, traj(grid))
}

plot_INLA = function(INLA_out, traj=NULL, xlim=NULL, ...)
{
  plot_BNPR(BNPR_out=INLA_out, traj, xlim, ...)
}

plot_INLA_inv = function(BNPR_out, traj=NULL, xlim=NULL, ...)
{
  mod = BNPR_out$result$summary.random$time
  
  grid = mod$"ID"
  if (is.null(xlim))
  {
    xlim=c(max(grid),0)
  }
  plot(grid,exp(mod$"0.5quant"),type="l",lwd=2.5,col="blue",log="y",
       xlab="Time (past to present)",ylab="Scaled Effective Pop. Size",
       xlim=xlim,ylim=c(min(exp(mod$"0.025quant"[grid > min(xlim) & grid < max(xlim)])),
                        max(exp(mod$"0.975quant"[grid > min(xlim) & grid < max(xlim)]))), ...)
  lines(grid,exp(mod$"0.975quant"),lwd=2.5,col="blue",lty=2)
  lines(grid,exp(mod$"0.025quant"),lwd=2.5,col="blue",lty=2)
  if (!is.null(traj))
    lines(grid, traj(grid))
}

plot_INLA_ii = function(BNPR_out, traj=NULL, xlim=NULL, ...)
{
  mod = BNPR_out$result$summary.random$ii
  
  grid = mod$"ID"
  if (is.null(xlim))
  {
    xlim=c(max(grid),0)
  }
  plot(grid,exp(-mod$"0.5quant"),type="l",lwd=2.5,col="blue",log="y",
       xlab="Time (past to present)",ylab="Scaled Effective Pop. Size",
       xlim=xlim,ylim=c(min(exp(-mod$"0.975quant"[grid > min(xlim) & grid < max(xlim)])),
                        max(exp(-mod$"0.025quant"[grid > min(xlim) & grid < max(xlim)]))), ...)
  lines(grid,exp(-mod$"0.975quant"),lwd=2.5,col="blue",lty=2)
  lines(grid,exp(-mod$"0.025quant"),lwd=2.5,col="blue",lty=2)
  if (!is.null(traj))
    lines(grid, traj(grid))
}
