#' Simulate from inhomogeneous, heterochronous coalescent
#' 
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled numeric vector of samples taken per sampling time.
#' @param traj function that returns effective population size at time t.
#' @param method which sampling method to use. "tt" invoke time-transformation
#'   method, "thin" invokes thinning method.
#' @param val_upper numeric used by time-transformation method to set a starting
#'   point for its dynamic numerical integration upper bound.
#' @param lower_bound numeric lower limit of \code{traj} function on its
#'   support.  Used only by thinning method.
#' @param ... additional arguments to be passed to \code{traj} function.
#'   
#' @return A list containing vectors of coalescent times \code{coal_times}, 
#'   intercoalescent times \code{intercoal_times}, and number of active lineages
#'   \code{lineages}, as well as passing along \code{samp_times} and
#'   \code{n_sampled}.
#' @export
#' 
#' @examples
#' coalsim(0:2, 3:1, unif_traj, lower_bound=10)
coalsim <- function(samp_times, n_sampled, traj, method="tt", val_upper=10, lower_bound=1, ...)
{
  if (method == "tt")
  {
    result = coalsim_tt(samp_times, n_sampled, traj, val_upper, ...)
  }
  else if (method == "thin")
  {
    result = coalsim_thin(samp_times, n_sampled, traj, lower_bound, ...)
  }
  else
  {
    stop("Argument method not recognized.")
  }
  
  return(result)
}

coalsim_tt <- function(samp_times, n_sampled, traj, val_upper=10, ...)
{
  traj_inv <- function(t) 1/traj(t, ...)
  hazard <- function(t, lins, start, target) .5*lins*(lins-1)*stats::integrate(traj_inv, start, start+t)$value - target
  
  coal_times = NULL
  lineages = NULL
  
  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]
  
  while (time <= max(samp_times) || active_lineages > 1)
  {
    if (active_lineages == 1)
    {
      curr <- curr + 1
      active_lineages <- active_lineages + n_sampled[curr]
      time <- samp_times[curr]
    }
    
    #time = time + stats::rexp(1, 0.5*active_lineages*(active_lineages-1)/lower_bound)
    target <- stats::rexp(1)
    y <- stats::uniroot(hazard, lins=active_lineages, start=time, target=target,
                 lower=0, upper=val_upper, extendInt = "upX")$root
    
    while(curr < length(samp_times) && time + y >= samp_times[curr+1])
    {
      target <- -hazard(t = samp_times[curr+1] - time, lins = active_lineages,
                        start = time, target = target)
      curr <- curr + 1
      active_lineages <- active_lineages + n_sampled[curr]
      time <- samp_times[curr]
      
      y <- stats::uniroot(hazard, lins=active_lineages, start=time, target=target,
                   lower=0, upper=val_upper, extendInt = "upX")$root
    }
    
    time <- time + y
    coal_times = c(coal_times, time)
    lineages = c(lineages, active_lineages)
    active_lineages = active_lineages - 1
  }
  
  return(list(coal_times = coal_times, lineages = lineages,
              intercoal_times = c(coal_times[1], diff(coal_times)),
              samp_times = samp_times, n_sampled = n_sampled))
}

coalsim_thin <- function(samp_times, n_sampled, traj, lower_bound, ...)
{
  coal_times = NULL
  lineages = NULL
  
  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]
  
  while (time <= max(samp_times) || active_lineages > 1)
  {
    if (active_lineages == 1)
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    
    time = time + stats::rexp(1, 0.5*active_lineages*(active_lineages-1)/lower_bound)
    
    if (curr < length(samp_times) && time >= samp_times[curr + 1])
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    else if (stats::runif(1) <= lower_bound/traj(time, ...))
    {
      coal_times = c(coal_times, time)
      lineages = c(lineages, active_lineages)
      active_lineages = active_lineages - 1
    }
  }
  
  return(list(coal_times = coal_times, lineages = lineages,
         intercoal_times = c(coal_times[1], diff(coal_times)),
         samp_times = samp_times, n_sampled = n_sampled))
}
