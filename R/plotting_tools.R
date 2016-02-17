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
                     lty=1, lwd=2, col="black", main="", log="y",
                     ylab="Effective Population Size",
                     xlab="Time", xmarline = 3, axlabs=NULL,
                     traj_lty=2, traj_lwd=2, traj_col=col,
                     newplot=TRUE, credible_region=TRUE,
                     heatmaps=TRUE, heatmap_labels=TRUE,
                     heatmap_labels_side="right", heatmap_width=7,
                     yscale = 1, ...)
{
  grid = BNPR_out$grid
  if (is.null(xlim))
  {
    xlim=c(max(grid), min(grid))
  }
  
  mask = BNPR_out$x >= min(xlim) & BNPR_out$x <= max(xlim)
  
  t = BNPR_out$x[mask]
  
  y = BNPR_out$effpop[mask] * yscale
  yhi = BNPR_out$effpop975[mask] * yscale
  ylo = BNPR_out$effpop025[mask] * yscale
  
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
      ylim=c(ymin/(yextra^1.35), ymax)
    }
    else
    {
      ylim=c(ymin, ymax)
    }
    
    if (is.null(axlabs))
    {
      plot(1,1,type="n",log=log,
           xlab=xlab, ylab=ylab, main=main,
           xlim=xlim, ylim=ylim, ...)
    }
    else
    {
      plot(1,1,type="n",log=log,
           xlab="", ylab=ylab, main=main,
           xlim=xlim, ylim=ylim, xaxt = "n", ...)
      axis(1, at=axlabs$x, labels = axlabs$labs, las=2)
      mtext(text = xlab, side = 1, line = xmarline)
    }
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

#' @export
plot_mrw = function(BNPR_outs, traj=NULL, xlim=NULL, nbreaks=40, ltys=1, lwds=2,
                    cols="black", xlab="Time", xmarline = 3, axlabs=NULL,
                    ylim=NULL, ymin_zero=FALSE, ylab="Mean Relative Width",
                    main="", heatmaps=TRUE, heatmap_labels=FALSE,
                    heatmap_labels_side="right", heatmap_width = 7, 
                    legends=NULL, legend_place="topleft",
                    legend_cex=1.0, bty="n", ...)
{
  n = length(BNPR_outs)
  ltys <- rep_len(x = ltys, length.out = n)
  lwds <- rep_len(x = lwds, length.out = n)
  cols <- rep_len(x = cols, length.out = n)
  
  x = BNPR_outs[[1]]$x
  if (is.null(xlim))
  {
    xlim=c(max(x),0)
  }
  
  mask = (x > min(xlim) & x < max(xlim))
  
  los  = list()
  his  = list()
  mids = list()
  mres = list()
  for (i in 1:length(BNPR_outs))
  {
    los[[i]]  = BNPR_outs[[i]]$effpop025[mask]
    his[[i]]  = BNPR_outs[[i]]$effpop975[mask]
    mids[[i]] = BNPR_outs[[i]]$effpop[mask]
    
    if (is.null(traj))
    {
      mres[[i]] = (his[[i]] - los[[i]])/mids[[i]]
    }
    else
    {
      mres[[i]] = (his[[i]] - los[[i]])/traj(x[mask])
    }
  }
  
  if (!is.null(ylim))
  {
    ymax = max(ylim)
    ymin = min(ylim)
  }
  else
  {
    ymax=max(unlist(mres))
    if (ymin_zero)
      ymin=0
    else
      ymin=min(unlist(mres))
  }
  
  if (heatmaps)
  {
    yspan=ymax-ymin
    yextra=yspan/10
    ylim=c(ymin-(yextra*1.35), ymax)
  }
  else
  {
    ylim = c(ymin, ymax)
  }
  
  if (is.null(axlabs))
  {
    plot(1,1,type="n",
         xlab=xlab,ylab=ylab,main=main,
         xlim=xlim,ylim=ylim)
  }
  else
  {
    plot(1,1,type="n", xaxt="n",
         xlab="",ylab=ylab,main=main,
         xlim=xlim,ylim=ylim, ...)
    axis(1, at=axlabs$x, labels = axlabs$labs, las=2, cex.lab=0.6)
    mtext(text = xlab, side = 1, line = xmarline)
  }
  
  for (i in 1:length(mres))
  {
    lines(x[mask], mres[[i]], lwd=lwds[i], col=cols[i], lty = ltys[i])
  }
  
  if (heatmaps)
  {
    samps = rep(BNPR_outs[[1]]$samp_times, BNPR_outs[[1]]$n_sampled)
    samps = samps[samps <= max(xlim) & samps >= min(xlim)]
    
    coals = BNPR_outs[[1]]$coal_times
    coals = coals[coals <= max(xlim) & coals >= min(xlim)]
    
    breaks = seq(min(xlim), max(xlim), length.out=nbreaks)
    h_samp = hist(samps, breaks=breaks, plot=FALSE)
    h_coal = hist(coals, breaks=breaks, plot=FALSE)
    
    hist2heat(h_samp, y=ymin-yextra*0.5, wd=heatmap_width)
    hist2heat(h_coal, y=ymin-yextra, wd=heatmap_width)
    
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
      
      text(x = lab_x, y = ymin-(yextra*0.20), labels = "Sampling events",
           adj = c(lab_adj, 0), cex = 0.7)
      text(x = lab_x, y = ymin-(yextra*1.25), labels = "Coalescent events",
           adj = c(lab_adj, 1), cex = 0.7)
    }
  }
  
  if (!is.null(legends))
    legend(legend_place, legends, lty=ltys, lwd = lwds, bty = bty, col=cols, cex=legend_cex)
}

#' @export
plot_seasonality = function(BNPR_out, zero_date, start = 0.0, years = NULL, 
                            period = 1.0, ylim = NULL, nbreaks = 40,
                            lty=1, lwd=2, col_years="gray", col_mean="black",
                            main="", axlabs = NULL, xlab = "Time", xmarline = 3,
                            ylab="Effective Population Size", log_y=TRUE,
                            heatmaps=TRUE, heatmap_labels=TRUE,
                            heatmap_labels_side="right", heatmap_width = 7,
                            legend = NULL, yscale = 1.0, ...)
{
  offset <- zero_date %% period - start
  if (log_y)
    log = "y"
  else
    log = ""
  
  t <- BNPR_out$x - offset
  y <- BNPR_out$effpop * yscale
  
  fun = approxfun(x = t, y = y, method = "linear")
  #yhi <- BNPR_out$effpop975[mask] * yscale
  #ylo <- BNPR_out$effpop025[mask] * yscale
  
  if (is.null(ylim))
  {
    ymax <- max(y)
    ymin <- min(y)
  }
  else
  {
    ymin <- min(ylim)
    ymax <- max(ylim)
  }
  
  if (heatmaps)
  {
    if (log_y)
    {
      yspan <- ymax/ymin
      yextra <- yspan^(1/10)
      ylim <- c(ymin/(yextra^1.35), ymax)
    }
    else
    {
      yspan <- ymax - ymin
      yextra <- yspan/10
      ylim <- c(ymin - yextra * 1.35, ymax)
    }
  }
  else
  {
    ylim=c(ymin, ymax)
  }
  
  if (is.null(axlabs))
  {
    plot(1, 1, type="n", log = log,
         xlab = xlab, ylab = ylab, main = main,
         xlim = c(period, 0), ylim = ylim, ...)
  }
  else
  {
    plot(1, 1, type="n", xaxt = "n", log = log,
         xlab = "", ylab = ylab, main = main,
         xlim = c(period, 0), ylim = ylim, ...)
    axis(1, at=axlabs$x, labels = axlabs$labs, las=2, cex.lab=1.0)
    mtext(text = xlab, side = 1, line = xmarline)
  }
  
  if (is.null(years))
  {
    years <- max(t)
  }
  
  current_t <- 0
  fun_sum <- rep(0, 101)
  fun_num <- 0
  while (current_t < years)
  {
    #period_mask <- (t >= current_t) & (t <= current_t + period)
    ts = seq(from = current_t, to = current_t + period, length.out = 101)
    
    lines(x = ts - current_t, y = fun(ts),
          lty = lty, lwd = lwd, col = col_years)
    
    if (fun_num > 0)
      fun_sum <- fun_sum + fun(ts)
    
    fun_num <- fun_num + 1
    
    current_t <- current_t + period
  }
  fun_mean <- fun_sum / fun_num
  lines(x = seq(0, 1, length.out = 101), y = fun_mean,
        lty = lty, lwd = lwd, col = col_mean)
  
  if (heatmaps)
  {
    samps = rep(BNPR_out$samp_times, BNPR_out$n_sampled) - offset
    samps = samps[samps <= years & samps >= 0]
    
    coals = BNPR_out$coal_times - offset
    coals = coals[coals <= years & coals >= 0]
    
    breaks = seq(0, period, length.out=nbreaks)
    h_samp = hist(samps %% period, breaks=breaks, plot=FALSE)
    h_coal = hist(coals %% period, breaks=breaks, plot=FALSE)
    
    if (log_y)
    {
      hist2heat(h_samp, y=ymin/yextra^0.5, wd=heatmap_width)
      hist2heat(h_coal, y=ymin/yextra, wd=heatmap_width)
    }
    else
    {
      hist2heat(h_samp, y=ymin - yextra * 0.5, wd=heatmap_width)
      hist2heat(h_coal, y=ymin - yextra, wd=heatmap_width)
    }
    
    #points(samps, rep(ymin/yextra^0.5, length(samps)), pch=3) 
    #points(coals, rep(ymin/yextra, length(coals)), pch=4)
    
    if (heatmap_labels)
    {
      if (heatmap_labels_side == "left")
      {
        lab_x = period
        lab_adj = 0
      }
      else if (heatmap_labels_side == "right")
      {
        lab_x = 0
        lab_adj = 1
      }
      else
      {
        warning('heatmap_labels_side not "left" or "right", defaulting to right')
        lab_x = 0
        lab_adj = 1
      }
      
      if (log_y)
      {
        text(x = lab_x, y = ymin/(yextra^0.20), labels = "Sampling events",
             adj = c(lab_adj, 0), cex = 0.7)
        text(x = lab_x, y = ymin/(yextra^1.25), labels = "Coalescent events",
             adj = c(lab_adj, 1), cex = 0.7)
      }
      else
      {
        text(x = lab_x, y = ymin - yextra * 0.20, labels = "Sampling events",
             adj = c(lab_adj, 0), cex = 0.7)
        text(x = lab_x, y = ymin - yextra * 1.25, labels = "Coalescent events",
             adj = c(lab_adj, 1), cex = 0.7)
      }
    }
  }
  
  if (!is.null(legend))
  {
    legend("topright", legend = legend, lty=1, lwd = 2, bty = "n", col=col_mean)
    
    #text(x=max(xlim), y=max(ylim), label=legend1, adj=c(0,1), cex=1.5)
  }
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
