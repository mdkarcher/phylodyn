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

#' @export
coalsim_tt <- function(samp_times, n_sampled, traj, val_upper=10, ...)
{
  traj_inv <- function(t) 1/traj(t, ...)
  hazard <- function(t, lins, start, target) .5*lins*(lins-1)*integrate(traj_inv, start, start+t)$value - target
  
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
    
    #time = time + rexp(1, 0.5*active_lineages*(active_lineages-1)/lower_bound)
    target <- rexp(1)
    y <- uniroot(hazard, lins=active_lineages, start=time, target=target,
                 lower=0, upper=val_upper, extendInt = "upX")$root
    
    while(curr < length(samp_times) && time + y >= samp_times[curr+1])
    {
      target <- -hazard(t = samp_times[curr+1] - time, lins = active_lineages,
                        start = time, target = target)
      curr <- curr + 1
      active_lineages <- active_lineages + n_sampled[curr]
      time <- samp_times[curr]
      
      y <- uniroot(hazard, lins=active_lineages, start=time, target=target,
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
    
    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1)/lower_bound)
    
    if (curr < length(samp_times) && time >= samp_times[curr + 1])
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    else if (runif(1) <= lower_bound/traj(time, ...))
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

#' Generate inhomogeneous coalescent with heterochronous samples.
#' 
#' @param sample a matrix with 2 columns. The first column contains the number 
#'   of samples collected at the time defined in the second column.
#' @param traj_inv a function returning the one over the effective population 
#'   size.
#' @param upper numeric. The maximum of the traj_inv function on the relevant 
#'   interval.
#' @param ... additional parameters to pass to the traj_inv function.
#'   
#' @return A list containing coalescent times \code{coal_times}, intercoalescent
#'   times \code{intercoal_times}, and number of active lineages between
#'   coalescent times \code{lineages}.
#' @export
#' 
#' @examples
#' coalgen_thinning_hetero(cbind(c(10,3), c(0,5)), unif_traj_inv, 0.01)
coalgen_thinning_hetero <- function(sample,traj_inv,upper,...)
{
  #'sample = is a matrix with 2 columns. The first column contains the number of samples collected at the time defined in the second column
  #'traj_inv = the effective population size function
  # this works for heterochronous sampling
  # assumes sample[1,1]>1
  samp_times = sample[,2]
  n_sampled = sample[,1]
  
  s=sample[1,2]
  b <- sample[1,1]
  n <- sum(sample[,1])-1
  m <- n
  nsample <- nrow(sample)
  sample <- rbind(sample,c(0,10*max(sample,2))) # Bug?
  out <- rep(0,n)
  branches <- rep(0,n)
  i <- 1
  while (i<(nsample))
  {
    #if (b==1)
    #{
    #  break
    #}
    if (b<2)
    {
      b <- b+sample[i+1,1]
      s <- sample[i+1,2]
      i <- i+1
    }
    E <- rexp(1,upper*b*(b-1)*.5)
    if (runif(1) <= traj_inv(E+s,...)/upper)
    {
      if ( (s+E)>sample[i+1,2])
      {
        b <- b+sample[i+1,1]
        s <- sample[i+1,2]
        i <- i+1
      }
      else
      {
        s <- s+E
        out[m-n+1] <- s
        branches[m-n+1] <- b
        n <- n-1
        b <- b-1
      }
    }
    else
    {
      s <- s+E
    }    
  }
  
  while (b>1)
  { 
    E <- rexp(1,upper*b*(b-1)*.5)
    if (runif(1)<=traj_inv(E+s,...)/upper)
    {
      s <- s+E
      out[m-n+1] <- s
      branches[m-n+1] <- b
      n <- n-1
      b <- b-1
    }
    else
    {
      s <- s+E
    }
  }
  
  return(list(intercoal_times=c(out[1],diff(out)), lineages=branches, 
              coal_times=out, samp_times = samp_times, n_sampled = n_sampled))   
}

coalgen_thinning_iso <- function(sample,traj_inv,upper=25,...)
{
  ###Need to add correction to "systematic" definition of upper bound
  s=sample[2]
  n <- sample[1]
  out <- rep(0,n-1)
  time <- 0
  j <- n
  while (j>1)
  {
    time <- time+rexp(1,upper*j*(j-1)*.5)
    if (runif(1)<=1/(traj_inv(time,...)*upper))
    {
      out[n-j+1] <- time
      j <- j-1
    }
  }
  return(list(intercoal_times=c(out[1],diff(out)),lineages=seq(n,2,-1), coal_times=out))
}

coalgen_transformation_hetero <- function(sample, trajectory,val_upper=10)
{
  #'sample = is a matrix with 2 columns. The first column contains the number of samples collected at the time defined in the second column
  #'trajectory = one over the effective population size function
  # this works for heterochronous sampling
  # assumes sample[1,1]>1
  s=sample[1,2]
  b <- sample[1,1]
  n <- sum(sample[,1])-1
  m <- n
  nsample <- nrow(sample)
  sample <- rbind(sample,c(0,10))
  out <- rep(0,n)
  branches <- rep(0,n)
  i <- 1
  while (i<(nsample+1))
  {
    #if (b==1)
    #{
    #  break
    #}
    if (b<2)
    {
      b <- b+sample[i+1,1]
      s <- sample[i+1,2]
      i <- i+1
    }
    x <- rexp(1)
    f <- function(bran,u,x,s) .5*bran*(bran-1)*integrate(trajectory, s, s+u)$value - x    
    y <- uniroot(f,bran=b,x=x,s=s,lower=0,upper=val_upper)$root
    while ( (s+y)>sample[i+1,2])
    {
      #     f <- function(bran,u,x,s) .5*bran*(bran-1)*integrate(trajectory, s, s+u)$value - x    
      #     y <- uniroot(f,bran=b,x=x,s=s,lower=0,upper=val_upper)$root
      x <- x-.5*b*(b-1)*integrate(trajectory,s,sample[i+1,2])$value  
      b <- b+sample[i+1,1]
      s <- sample[i+1,2]
      i <- i+1
      f <- function(bran,u,x,s) .5*bran*(bran-1)*integrate(trajectory, s, s+u)$value - x    
      y <- uniroot(f,bran=b,x=x,s=s,lower=0,upper=val_upper)$root
      if (i==nsample)
      {
        sample[nsample+1,2] <- 10*(s+y)
      }
    } 
    
    s <- s+y
    out[m-n+1] <- s
    branches[m-n+1] <- b
    n <- n-1
    b <- b-1
    if (i==nsample)
    {
      sample[nsample+1,2] <- 10*(s+y)
    }
  }
  
  return(list(branches=c(out[1],diff(out)),lineages=branches))   
}