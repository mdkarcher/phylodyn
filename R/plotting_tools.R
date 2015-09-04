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

#' @export
plot_BNPR = function(BNPR_out, traj=NULL, xlim=NULL, ylim=NULL, nbreaks=40,
                     lty=1, lwd=2, col="black", main="", xlab="Time",
                     ylab="Effective Population Size", log="y",
                     traj_lty=2, traj_lwd=2, traj_col=col,
                     newplot=TRUE, credible_region=TRUE,
                     heatmaps=TRUE, heatmap_labels=TRUE,
                     heatmap_labels_side="right",
                     heatmap_width=10,...)
{
  grid = BNPR_out$grid
  if (is.null(xlim))
  {
    xlim=c(max(grid), min(grid))
  }
  
  mask = BNPR_out$x >= min(xlim) & BNPR_out$x <= max(xlim)
  
  t = BNPR_out$x[mask]
  
  y = BNPR_out$effpop[mask]
  yhi = BNPR_out$effpop975[mask]
  ylo = BNPR_out$effpop025[mask]
  
  if (newplot)
  {
    if (is.null(ylim))
    {
      ymax=max(yhi)
      ymin=min(ylo)
    }
    else
    {
      ymin=min(ylim)
      ymax=max(ylim)
    }
    
    if (heatmaps)
    {
      yspan=ymax/ymin
      yextra=yspan^(1/10)
      ylim=c(ymin/(yextra^1.25), ymax)
    }
    else
    {
      ylim=c(ymin, ymax)
    }
    
    plot(1,1,type="n",log=log,
         xlab=xlab, ylab=ylab, main=main,
         xlim=xlim, ylim=ylim, ...)
  }
  
  if (credible_region)
  {
    shade_band(x = t, ylo = ylo, yhi = yhi, col="lightgray")
  }
  
  if (!is.null(traj))
  {
    lines(t, traj(t), lwd=traj_lwd, lty=traj_lty, col=traj_col)
  }
  
  if (newplot)
  {
    if (heatmaps)
    {
      samps = rep(BNPR_out$samp_times, BNPR_out$n_sampled)
      samps = samps[samps <= max(xlim) & samps >= min(xlim)]
      
      coals = BNPR_out$coal_times
      coals = coals[coals <= max(xlim) & coals >= min(xlim)]
      
      breaks = seq(min(xlim), max(xlim), length.out=nbreaks)
      h_samp = hist(samps, breaks=breaks, plot=FALSE)
      h_coal = hist(coals, breaks=breaks, plot=FALSE)
      
      hist2heat(h_samp, y=ymin/yextra^0.5, wd=heatmap_width)
      hist2heat(h_coal, y=ymin/yextra, wd=heatmap_width)
      
      #points(samps, rep(ymin/yextra^0.5, length(samps)), pch=3) 
      #points(coals, rep(ymin/yextra, length(coals)), pch=4)
      
      if (heatmap_labels)
      {
        if (heatmap_labels_side == "left")
        {
          lab_x = max(xlim)
          lab_adj = 0
        }
        else if (heatmap_labels_side == "right")
        {
          lab_x = min(xlim)
          lab_adj = 1
        }
        else
        {
          warning('heatmap_labels_side not "left" or "right", defaulting to right')
          lab_x = min(xlim)
          lab_adj = 1
        }
        
        text(x = lab_x, y = ymin/(yextra^0.20), labels = "Sampling events",
             adj = c(lab_adj, 0), cex = 0.7)
        text(x = lab_x, y = ymin/(yextra^1.25), labels = "Coalescent events",
             adj = c(lab_adj, 1), cex = 0.7)
      }
    }
  }
  
  lines(t, y, lwd=lwd, col=col, lty=lty)
}

plot_BNPR_old = function(BNPR_out, traj=NULL, xlim=NULL, ...)
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
