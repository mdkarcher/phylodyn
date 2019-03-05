shade_band = function(x, ylo, yhi, xlim=NULL, col="gray")
{
  if (is.null(xlim))
    xlim = c(0, Inf)
  mask = x >= min(xlim) & x <= max(xlim)
  
  x = x[mask]
  ylo = ylo[mask]
  yhi = yhi[mask]
  
  graphics::polygon(c(x, rev(x)), c(yhi, rev(ylo)), col=col, border=NA)
}

#' Plot a BNPR output
#' 
#' @param BNPR_out output of BNPR or BNPR_PS.
#' @param traj function summarizing the true effective population size 
#'   trajectory.
#' @param xlim numeric x-axis interval.
#' @param ylim numeric y-axis interval.
#' @param nbreaks integer number of bins for sampling heatmap.
#' @param lty line type for estimated trajectory.
#' @param lwd line width for estimated trajectory.
#' @param col color for estimated trajectory.
#' @param main character main plot title.
#' @param log character which axes to plot log-scale. Defaults to "y".
#' @param ylab character y-axis label.
#' @param xlab character x-axis label.
#' @param xmarline numeric if not using default x-axis labels, how far to put 
#'   the labels from the axis.
#' @param axlabs character vector x-axis labels.
#' @param traj_lty,traj_lwd,traj_col line type, line width, and line color for 
#'   the true trajectory.
#' @param newplot boolean whether to create a new plot or superimpose over a 
#'   previously open plot.
#' @param credible_region logical whether to display pointwise credible region.
#' @param heatmaps boolean whether to display sampling and coalescent heatmaps.
#' @param heatmap_labels boolean whether to display labels on heatmaps.
#' @param heatmap_labels_side string which side of plot to display heatmaps.
#' @param heatmap_labels_cex numeric scaling factor for heatmap labels.
#' @param heatmap_width numeric how wide heatmaps should be.
#' @param yscale numeric scaling applied to all effective population
#'   calculations.
#' @param ... additional arguments to be passed onto plot().
#'   
#' @export
plot_BNPR = function(BNPR_out, traj=NULL, xlim=NULL, ylim=NULL, nbreaks=40,
                     lty=1, lwd=2, col="black", main="", log="y",
                     ylab="Effective Population Size",
                     xlab="Time", xmarline = 3, axlabs=NULL,
                     traj_lty=2, traj_lwd=2, traj_col=col,
                     newplot=TRUE, credible_region=TRUE,
                     heatmaps=TRUE, heatmap_labels=TRUE,
                     heatmap_labels_side="right", heatmap_labels_cex = 0.7,
                     heatmap_width=7, yscale = 1, ...)
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
      graphics::plot(1,1,type="n",log=log,
                     xlab=xlab, ylab=ylab, main=main,
                     xlim=xlim, ylim=ylim, ...)
    }
    else
    {
      graphics::plot(1,1,type="n",log=log,
                     xlab="", ylab=ylab, main=main,
                     xlim=xlim, ylim=ylim, xaxt = "n", ...)
      graphics::axis(1, at=axlabs$x, labels = axlabs$labs, las=2)
      graphics::mtext(text = xlab, side = 1, line = xmarline)
    }
  }
  
  if (credible_region)
  {
    shade_band(x = t, ylo = ylo, yhi = yhi, col="lightgray")
  }
  
  if (!is.null(traj))
  {
    graphics::lines(t, traj(t), lwd=traj_lwd, lty=traj_lty, col=traj_col)
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
      h_samp = graphics::hist(samps, breaks=breaks, plot=FALSE)
      h_coal = graphics::hist(coals, breaks=breaks, plot=FALSE)
      
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
        
        graphics::text(x = lab_x, y = ymin/(yextra^0.20), labels = "Sampling events",
             adj = c(lab_adj, 0), cex = heatmap_labels_cex)
        graphics::text(x = lab_x, y = ymin/(yextra^1.25), labels = "Coalescent events",
             adj = c(lab_adj, 1), cex = heatmap_labels_cex)
      }
    }
  }
  
  graphics::lines(t, y, lwd=lwd, col=col, lty=lty)
}

#' Plot Mean Relative Widths
#' 
#' @param BNPR_outs one or more outputs of BNPR or BNPR_PS.
#' @param traj function summarizing the true effective population size 
#'   trajectory.
#' @param xlim numeric x-axis interval.
#' @param ylim numeric y-axis interval.
#' @param nbreaks integer number of bins for sampling heatmap.
#' @param ltys line types for estimated trajectories.
#' @param lwds line widths for estimated trajectories.
#' @param cols colors for estimated trajectories.
#' @param xlab character x-axis label.
#' @param xmarline numeric if not using default x-axis labels, how far to put
#'   the labels from the axis.
#' @param axlabs character vector x-axis labels.
#' @param ymin_zero logical 
#' @param ylab character y-axis label.
#' @param main character main plot title.
#' @param heatmaps boolean whether to display sampling and coalescent heatmaps.
#' @param heatmap_labels boolean whether to display labels on heatmaps.
#' @param heatmap_labels_side string which side of plot to display heatmaps.
#' @param heatmap_width numeric how wide heatmaps should be.
#' @param legends character legend texts.
#' @param legend_place character location of legend pane. See legend().
#' @param legend_cex numeric expansion factor for legend pane.
#' @param bty integer box type. See legend().
#' @param ... additional arguments to be passed onto plot().
#'
#' @export
plot_mrw = function(BNPR_outs, traj=NULL, xlim=NULL, ylim=NULL, nbreaks=40,
                    ltys=1, lwds=2, cols="black", xlab="Time", xmarline = 3,
                    axlabs=NULL, ymin_zero=FALSE, ylab="Mean Relative Width",
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
    graphics::plot(1,1,type="n",
                   xlab=xlab,ylab=ylab,main=main,
                   xlim=xlim,ylim=ylim)
  }
  else
  {
    graphics::plot(1,1,type="n", xaxt="n",
                   xlab="",ylab=ylab,main=main,
                   xlim=xlim,ylim=ylim, ...)
    graphics::axis(1, at=axlabs$x, labels = axlabs$labs, las=2, cex.lab=0.6)
    graphics::mtext(text = xlab, side = 1, line = xmarline)
  }
  
  for (i in 1:length(mres))
  {
    graphics::lines(x[mask], mres[[i]], lwd=lwds[i], col=cols[i], lty = ltys[i])
  }
  
  if (heatmaps)
  {
    samps = rep(BNPR_outs[[1]]$samp_times, BNPR_outs[[1]]$n_sampled)
    samps = samps[samps <= max(xlim) & samps >= min(xlim)]
    
    coals = BNPR_outs[[1]]$coal_times
    coals = coals[coals <= max(xlim) & coals >= min(xlim)]
    
    breaks = seq(min(xlim), max(xlim), length.out=nbreaks)
    h_samp = graphics::hist(samps, breaks=breaks, plot=FALSE)
    h_coal = graphics::hist(coals, breaks=breaks, plot=FALSE)
    
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
      
      graphics::text(x = lab_x, y = ymin-(yextra*0.20), labels = "Sampling events",
           adj = c(lab_adj, 0), cex = 0.7)
      graphics::text(x = lab_x, y = ymin-(yextra*1.25), labels = "Coalescent events",
           adj = c(lab_adj, 1), cex = 0.7)
    }
  }
  
  if (!is.null(legends))
    graphics::legend(legend_place, legends, lty=ltys, lwd = lwds, bty = bty, col=cols, cex=legend_cex)
}

#' Plot Seasonality
#' 
#' @param BNPR_out output of BNPR or BNPR_PS.
#' @param zero_date numeric which time to present as 0 on graph.
#' @param start numeric when to start displaying data.
#' @param years numeric how many years to count back in time.
#' @param period numeric how long of a period is a "year".
#' @param ylim numeric y-axis interval.
#' @param nbreaks integer number of bins for sampling heatmap.
#' @param lty numeric line type for the estimated trajectory paths.
#' @param lwd numeric line width for the estimated trajectory paths.
#' @param col_years color of year paths.
#' @param col_mean color of mean path.
#' @param main character main plot title.
#' @param axlabs character vector x-axis labels.
#' @param xlab character x-axis label.
#' @param xmarline numeric if not using default x-axis labels, how far to put
#'   the labels from the axis.
#' @param ylab character y-axis label.
#' @param log_y logical whether to log-scale the y-axis.
#' @param heatmaps boolean whether to display sampling and coalescent heatmaps.
#' @param heatmap_labels boolean whether to display labels on heatmaps.
#' @param heatmap_labels_side string which side of plot to display heatmaps.
#' @param heatmap_width numeric how wide heatmaps should be.
#' @param legend character legend text. The defauly, NULL, disables.
#' @param yscale numeric scaling applied to all effective population
#'   calculations.
#' @param ... additional arguments to be passed onto plot().
#'
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
  
  fun = stats::approxfun(x = t, y = y, method = "linear")
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
    graphics::plot(1, 1, type="n", log = log,
                   xlab = xlab, ylab = ylab, main = main,
                   xlim = c(period, 0), ylim = ylim, ...)
  }
  else
  {
    graphics::plot(1, 1, type="n", xaxt = "n", log = log,
                   xlab = "", ylab = ylab, main = main,
                   xlim = c(period, 0), ylim = ylim, ...)
    graphics::axis(1, at=axlabs$x, labels = axlabs$labs, las=2, cex.lab=1.0)
    graphics::mtext(text = xlab, side = 1, line = xmarline)
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
    
    graphics::lines(x = ts - current_t, y = fun(ts),
          lty = lty, lwd = lwd, col = col_years)
    
    if (fun_num > 0)
      fun_sum <- fun_sum + fun(ts)
    
    fun_num <- fun_num + 1
    
    current_t <- current_t + period
  }
  fun_mean <- fun_sum / fun_num
  graphics::lines(x = seq(0, 1, length.out = 101), y = fun_mean,
        lty = lty, lwd = lwd, col = col_mean)
  
  if (heatmaps)
  {
    samps = rep(BNPR_out$samp_times, BNPR_out$n_sampled) - offset
    samps = samps[samps <= years & samps >= 0]
    
    coals = BNPR_out$coal_times - offset
    coals = coals[coals <= years & coals >= 0]
    
    breaks = seq(0, period, length.out=nbreaks)
    h_samp = graphics::hist(samps %% period, breaks=breaks, plot=FALSE)
    h_coal = graphics::hist(coals %% period, breaks=breaks, plot=FALSE)
    
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
        graphics::text(x = lab_x, y = ymin/(yextra^0.20), labels = "Sampling events",
             adj = c(lab_adj, 0), cex = 0.7)
        graphics::text(x = lab_x, y = ymin/(yextra^1.25), labels = "Coalescent events",
             adj = c(lab_adj, 1), cex = 0.7)
      }
      else
      {
        graphics::text(x = lab_x, y = ymin - yextra * 0.20, labels = "Sampling events",
             adj = c(lab_adj, 0), cex = 0.7)
        graphics::text(x = lab_x, y = ymin - yextra * 1.25, labels = "Coalescent events",
             adj = c(lab_adj, 1), cex = 0.7)
      }
    }
  }
  
  if (!is.null(legend))
  {
    graphics::legend("topright", legend = legend, lty=1, lwd = 2, bty = "n", col=col_mean)
    
    #graphics::text(x=max(xlim), y=max(ylim), label=legend1, adj=c(0,1), cex=1.5)
  }
}

plot_MCMC = function(MCMC_out, traj=NULL, xlim=NULL, ylim=NULL, nbreaks=40,
                     lty=1, lwd=2, col="black", main="", log="y",
                     ylab="Effective Population Size",
                     xlab="Time", xmarline = 3, axlabs=NULL,
                     traj_lty=2, traj_lwd=2, traj_col=col,
                     newplot=TRUE, credible_region=TRUE,
                     heatmaps=TRUE, heatmap_labels=TRUE,
                     heatmap_labels_side="right", heatmap_width=7,
                     yscale = 1, ...)
{
  grid = MCMC_out$grid
  if (is.null(xlim))
  {
    xlim=c(max(grid), min(grid))
  }
  
  mask = MCMC_out$x >= min(xlim) & MCMC_out$x <= max(xlim)
  
  t = MCMC_out$x[mask]
  
  y = MCMC_out$med_fun(t) * yscale
  yhi = MCMC_out$hi_fun(t) * yscale
  ylo = MCMC_out$low_fun(t) * yscale
  
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
      graphics::plot(1,1,type="n",log=log,
                     xlab=xlab, ylab=ylab, main=main,
                     xlim=xlim, ylim=ylim, ...)
    }
    else
    {
      graphics::plot(1,1,type="n",log=log,
                     xlab="", ylab=ylab, main=main,
                     xlim=xlim, ylim=ylim, xaxt = "n", ...)
      graphics::axis(1, at=axlabs$x, labels = axlabs$labs, las=2)
      graphics::mtext(text = xlab, side = 1, line = xmarline)
    }
  }
  
  if (credible_region)
  {
    shade_band(x = t, ylo = ylo, yhi = yhi, col="lightgray")
  }
  
  if (!is.null(traj))
  {
    graphics::lines(t, traj(t), lwd=traj_lwd, lty=traj_lty, col=traj_col)
  }
  
  if (newplot)
  {
    if (heatmaps)
    {
      samps = rep(MCMC_out$samp_times, MCMC_out$n_sampled)
      samps = samps[samps <= max(xlim) & samps >= min(xlim)]
      
      coals = MCMC_out$coal_times
      coals = coals[coals <= max(xlim) & coals >= min(xlim)]
      
      breaks = seq(min(xlim), max(xlim), length.out=nbreaks)
      h_samp = graphics::hist(samps, breaks=breaks, plot=FALSE)
      h_coal = graphics::hist(coals, breaks=breaks, plot=FALSE)
      
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
        
        graphics::text(x = lab_x, y = ymin/(yextra^0.20), labels = "Sampling events",
             adj = c(lab_adj, 0), cex = 0.7)
        graphics::text(x = lab_x, y = ymin/(yextra^1.25), labels = "Coalescent events",
             adj = c(lab_adj, 1), cex = 0.7)
      }
    }
  }
  
  graphics::lines(t, y, lwd=lwd, col=col, lty=lty)
}

plot_BNPR_old = function(BNPR_out, traj=NULL, xlim=NULL, ...)
{
  mod = BNPR_out$result$summary.random$time
  
  grid = mod$"ID"
  if (is.null(xlim))
  {
    xlim=c(max(grid),0)
  }
  graphics::plot(grid,exp(-mod$"0.5quant"),type="l",lwd=2.5,col="blue",log="y",
                 xlab="Time (past to present)",ylab="Scaled Effective Pop. Size",
                 xlim=xlim,ylim=c(min(exp(-mod$"0.975quant"[grid > min(xlim) & grid < max(xlim)])),
                                  max(exp(-mod$"0.025quant"[grid > min(xlim) & grid < max(xlim)]))), ...)
  graphics::lines(grid,exp(-mod$"0.975quant"),lwd=2.5,col="blue",lty=2)
  graphics::lines(grid,exp(-mod$"0.025quant"),lwd=2.5,col="blue",lty=2)
  if (!is.null(traj))
    graphics::lines(grid, traj(grid))
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
  graphics::plot(grid,exp(mod$"0.5quant"),type="l",lwd=2.5,col="blue",log="y",
                 xlab="Time (past to present)",ylab="Scaled Effective Pop. Size",
                 xlim=xlim,ylim=c(min(exp(mod$"0.025quant"[grid > min(xlim) & grid < max(xlim)])),
                                  max(exp(mod$"0.975quant"[grid > min(xlim) & grid < max(xlim)]))), ...)
  graphics::lines(grid,exp(mod$"0.975quant"),lwd=2.5,col="blue",lty=2)
  graphics::lines(grid,exp(mod$"0.025quant"),lwd=2.5,col="blue",lty=2)
  if (!is.null(traj))
    graphics::lines(grid, traj(grid))
}

plot_INLA_ii = function(BNPR_out, traj=NULL, xlim=NULL, ...)
{
  mod = BNPR_out$result$summary.random$ii
  
  grid = mod$"ID"
  if (is.null(xlim))
  {
    xlim=c(max(grid),0)
  }
  graphics::plot(grid,exp(-mod$"0.5quant"),type="l",lwd=2.5,col="blue",log="y",
                 xlab="Time (past to present)",ylab="Scaled Effective Pop. Size",
                 xlim=xlim,ylim=c(min(exp(-mod$"0.975quant"[grid > min(xlim) & grid < max(xlim)])),
                                  max(exp(-mod$"0.025quant"[grid > min(xlim) & grid < max(xlim)]))), ...)
  graphics::lines(grid,exp(-mod$"0.975quant"),lwd=2.5,col="blue",lty=2)
  graphics::lines(grid,exp(-mod$"0.025quant"),lwd=2.5,col="blue",lty=2)
  if (!is.null(traj))
    graphics::lines(grid, traj(grid))
}


#' Plot trajectory estimate and credible intervals
#'
#' @param midpts numeric vector of x-values at which to plot estimate and credible interval bounds.
#' @param med numeric vector of estimated trajectory values.
#' @param up numeric vector of trajectory pointwise credible interval upper bounds.
#' @param low numeric vector of trajectory pointwise credible interval lower bounds.
#' @param cutoff numeric maximum time (in past) plotted.
#' @param ylim numeric y-axis interval.
#' @param xlab character x-axis label.
#' @param ylab character y-axis label.
#' @param axlabs character vector x-axis labels.
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled integer vector of number of samples taken at each element of samp_times.
#' @param nbreaks integer number of bins in the sampling time histogram heatmap.
#' @param med_col color of trajectory estimate line.
#' @param cred_col color of trajectory pointwise credible interval region.
#' @param heatmap_width numeric vertical width of sampling time histogram heatmap.
#' @param heatmap_legend_cex numeric plotting expansion for heatmap legend.
#' @param ... additional parameters to be passed to plot().
#'
#' @return
#' @export
#'
#' @examples
plot_BEAST = function(midpts, med, up, low, cutoff = NULL, ylim = NULL,
                      xlab = "Time", ylab = "Effective Population Size",axlabs=NULL,
                      samp_times = NULL, n_sampled = NULL, nbreaks = 40, 
                      med_col = rgb(0.330, 0.484, 0.828), cred_col = rgb(0.330, 0.484, 0.828, 0.4),
                      heatmap_width = 7, heatmap_legend_cex = 1.0, ...) {
  #trev_blue = rgb(0.330, 0.484, 0.828)
  #trev_yell = rgb(0.829, 0.680, 0.306)
  
  if (is.null(cutoff)) {
    cutoff = max(midpts)
  }
  
  mask = midpts <= cutoff
  
  if (is.null(ylim)) {
    ylim = c(min(low[mask]), max(up[mask]))
  } else {
    ylim = c(min(ylim), max(ylim))
  }
  
  if (!is.null(samp_times)) {
    yspan=max(ylim)/min(ylim)
    yextra=yspan^(0.07)
    ylim[1] = ylim[1] / yextra
  }
  
  xaxt = "s"
  if (!is.null(axlabs))
  {
    xaxt = "n"
  }
  
  plot(x = midpts[mask], y = med[mask], 
       xlab = "", ylab = "", xlim = c(cutoff, 0), ylim = ylim, 
       log = "y", xaxt = xaxt, type = 'n', ...)
  
  if (!is.null(axlabs))
  {
    graphics::axis(1, at=axlabs$x, labels = axlabs$labs, las=0)
    #graphics::mtext(text = xlab, side = 1, line = 3)
  }
  
  if (!is.null(samp_times)) {
    if (is.null(n_sampled)) {
      n_sampled = rep(1, length(samp_times))
    }
    samps = rep(samp_times, n_sampled)
    samps = samps[samps <= cutoff]
    
    breaks = seq(0, cutoff, length.out=nbreaks)
    h_samp = graphics::hist(samps, breaks=breaks, plot=FALSE)
    
    hist2heat(h_samp, y=ylim[1], wd=heatmap_width)
    
    lab_x = 0
    lab_adj = 1
    
    graphics::text(x = lab_x, y = ylim[1]*yextra, labels = "Sampling events",
                   adj = c(1, 1), cex = heatmap_legend_cex)
  }
  
  phylodyn:::shade_band(x = midpts[mask], ylo = low[mask], yhi = up[mask],
                        col = cred_col)
  
  lines(x = midpts[mask], y = med[mask], lwd = 2, col = med_col)
  
  title(ylab = ylab, line = 2.5)
  title(xlab = xlab, line = 2.5)
}

#' Plot seasonal overlays
#'
#' @param BEAST_out list containing `midpts`, `medianEffPop`, `q975EffPop`, and `q025EffPop`.
#' @param zero_date numeric date (in decimal years) representing `t=0` in `BEAST_out`
#' @param start numeric value between 0 and `period` controlling the right side of the plot.
#' @param years numeric how many years or periods to overlay.
#' @param period numeric representing how long one year or other seasonal period is in the units of `BEAST_out`.
#' @param ylim numeric y-axis interval.
#' @param nbreaks integer number of bins in the sampling time histogram heatmap.
#' @param lty numeric line type.
#' @param lwd numeric line width.
#' @param col_years color of individual trajectory
#' @param col_median color of the median estimate.
#' @param main character main figure title.
#' @param axlabs 
#' @param xlab 
#' @param xmarline 
#' @param ylab 
#' @param log_y 
#' @param samp_times 
#' @param n_sampled 
#' @param heatmaps 
#' @param heatmap_labels 
#' @param heatmap_labels_side 
#' @param heatmap_width numeric vertical width of sampling time histogram heatmap.
#' @param legend 
#' @param yscale 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_seasonality = function(BEAST_out, zero_date, start = 0.0, years = NULL, 
                            period = 1.0, ylim = NULL, nbreaks = 40,
                            lty=1, lwd=2, col_years="gray", col_median="black",
                            main="", axlabs = NULL, xlab = "Time", xmarline = 1,
                            ylab="Effective Population Size", log_y=TRUE,
                            samp_times = NULL, n_sampled = NULL,
                            heatmaps=TRUE, heatmap_labels=TRUE,
                            heatmap_labels_side="right", heatmap_width = 7,
                            legend = NULL, yscale = 1.0, ...)
{
  offset <- zero_date %% period - start
  if (log_y)
    log = "y"
  else
    log = ""
  
  if (is.null(samp_times))
    heatmaps = FALSE
  
  t <- BEAST_out$midpts - offset
  y <- BEAST_out$medianEffPop * yscale
  
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
      yextra <- yspan^(0.05)
      ylim <- c(ymin/(yextra^1), ymax)
    }
    else
    {
      yspan <- ymax - ymin
      yextra <- yspan*0.05
      ylim <- c(ymin - yextra * 1, ymax)
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
    axis(1, at=axlabs$x, labels = axlabs$labs, las=2)
    mtext(text = xlab, side = 1, line = xmarline, cex = 1.3)
  }
  
  if (is.null(years))
  {
    years <- max(t)
  }
  
  current_t <- 0
  fun_mat <- matrix(0, nrow = 101, ncol = 0)
  fun_col <- 0
  while (current_t < years)
  {
    #period_mask <- (t >= current_t) & (t <= current_t + period)
    ts = seq(from = current_t, to = current_t + period, length.out = 101)
    
    lines(x = ts - current_t, y = fun(ts),
          lty = lty, lwd = lwd, col = col_years)
    
    if (fun_col > 0)
      fun_mat <- cbind(fun_mat, fun(ts))
    
    fun_col <- fun_col + 1
    
    current_t <- current_t + period
  }
  fun_median <- apply(X = fun_mat, MARGIN = 1, FUN = median)
  lines(x = seq(0, 1, length.out = 101), y = fun_median,
        lty = lty, lwd = lwd, col = col_median)
  
  if (heatmaps)
  {
    samps = rep(samp_times, n_sampled) - offset
    samps = samps[samps <= years & samps >= 0]
    
    #coals = BEAST_out$coal_times - offset
    #coals = coals[coals <= years & coals >= 0]
    
    breaks = seq(0, period, length.out=nbreaks)
    h_samp = hist(samps %% period, breaks=breaks, plot=FALSE)
    #h_coal = hist(coals %% period, breaks=breaks, plot=FALSE)
    
    if (log_y)
    {
      hist2heat(h_samp, y=ymin/yextra^1, wd=heatmap_width)
      #hist2heat(h_coal, y=ymin/yextra, wd=heatmap_width)
    }
    else
    {
      hist2heat(h_samp, y=ymin - yextra * 1, wd=heatmap_width)
      #hist2heat(h_coal, y=ymin - yextra, wd=heatmap_width)
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
        text(x = lab_x, y = ymin/(yextra^0.0), labels = "Sampling events",
             adj = c(lab_adj, 0), cex = 1)
        #text(x = lab_x, y = ymin/(yextra^1.25), labels = "Coalescent events",
        #     adj = c(lab_adj, 1), cex = 0.7)
      }
      else
      {
        text(x = lab_x, y = ymin - yextra * 0.0, labels = "Sampling events",
             adj = c(lab_adj, 0), cex = 1)
        #text(x = lab_x, y = ymin - yextra * 1.25, labels = "Coalescent events",
        #     adj = c(lab_adj, 1), cex = 0.7)
      }
    }
  }
  
  if (!is.null(legend))
  {
    legend("topright", legend = legend, lty=1, lwd = 2, bty = "n", col=col_median)
    
    #text(x=max(xlim), y=max(ylim), label=legend1, adj=c(0,1), cex=1.5)
  }
}
