#' Raise a trajectory function to a power
#' 
#' @param t numeric vector of times at which to evaluate the transformed
#'   trajectory function.
#' @param traj function. The trajectory to transform.
#' @param beta numeric. The power to raise the trajectory function to.
#' @param ... additional arguments to pass to trajectory function.
#'   
#' @return A vector containing \eqn{f(t)^\beta}.
#' @export
#' 
#' @examples
#' traj_beta(1:3, unif_traj, 2, level=2)
traj_beta = function(t, traj, beta, ...)
{
  return(traj(t, ...)^beta)
}

#' Uniform trajectory.
#' 
#' @param t numeric vector of times at which to evaluate the trajectory 
#'   function.
#' @param level numeric. The value to return at all values \code{t}.
#'   
#' @return A vector containing a value of \code{level} for every element in
#'   \code{t}.
#' @export
#' 
#' @examples
#' unif_traj(0:10, level=5)
unif_traj = function(t, level=100)
{
  n = length(t)
  return(rep(level,n))
}

#' @export
unif_traj_inv = function(t, level=100)
{
  n = length(t)
  return(rep(1/level,n))
}

#' Exponential trajectory
#'
#' @param t numeric vector of times at which to evaluate the trajectory 
#'   function.
#' @param scale value to return at \code{t=0}
#'
#' @return A vector containing a value of \code{level} for every element in
#'   \code{t}.
#' @export
#'
#' @examples
#' unif_traj(0:10, scale = 100)
exp_traj = function(t, scale=1000)
{
  return(scale * exp(-t))
}

#' @export
exp_traj_inv = function(t, scale=1000)
{
  return(1/(scale * exp(-t)))
}

#' @export
boombust_traj = function(t, bust=1, scale=1000)
{
  #bust = 1
  result = rep(0, length(t))
  result[t <= bust] = scale*exp(t[t <= bust]-bust)
  result[t >  bust] = scale*exp(bust-t[t >  bust])
  return(result)
}

#' @export
boombust_traj_inv = function(t, bust=1, scale=1000)
{
  return(1/boombust_traj(t, bust, scale))
}

#' @export
cyclic_traj = function(t)
{
  result = rep(0, length(t))
  result[(t %% 10) <= 5] = 200*exp(-(t[(t %% 10) <= 5] %% 10) / 2)
  result[(t %% 10) >  5] = 200*exp((t[(t %% 10) >  5] %% 10) / 2 - 5)
  return(result)
}

#' @export
cyclic_traj_inv = function(t)
{
  return(1/cyclic_traj(t))
}

#' @export
steep_cyc_traj = function(t)
{
  result = rep(0, length(t))
  result[(t %% 10) <= 5] = 20 + 1980*exp(-(t[(t %% 10) <= 5] %% 10) * 5)
  result[(t %% 10) >  5] = 20 + 1980*exp(((t[(t %% 10) >  5] %% 10) - 10) * 5)
  return(result)
}

#' @export
steep_cyc_traj_inv = function(t)
{
  return(1/steep_cyc_traj(t))
}

#' @export
sloped_traj = function(t)
{
  result = rep(0, length(t))
  result[(t %% 10) <= 5] = 100 + 900*exp(-(t[(t %% 10) <= 5] %% 10) * 1)
  result[(t %% 10) >  5] = 100 + 900*exp(((t[(t %% 10) >  5] %% 10) - 10) * 1)
  return(result)
}

#' @export
sloped_traj_inv = function(t)
{
  return(1/sloped_traj(t))
}

#' @export
mesa_traj = function(t, a=2, b=3)
{
  result = rep(0, length(t))
  result[(t %% 10) <= a | (t %% 10) > 10-a] = 1000
  result[(t %% 10) >  a & (t %% 10) <= 5] = 100 + 900*exp((a-((t[(t %% 10) >  a & (t %% 10) <= 5] %% 10))) * b)
  result[(t %% 10) >  5 & (t %% 10) <= 10-a] = 100 + 900*exp(((t[(t %% 10) >  5 & (t %% 10) <= 10-a] %% 10) - 10 + a) * b)
  return(result)
}

#' @export
mesa_traj_inv = function(t, a=2, b=3)
{
  return(1/mesa_traj(t, a, b))
}

#' Seasonal logistic trajectory
#'
#' @param t numeric vector of times at which to evaluate the trajectory 
#'   function.
#' @param offset numeric. By default \code{t=0} is a local maximum. Offset is
#'   added to \code{t} to allow \code{t=0} to represent a different point in the
#'   cycle.
#' @param a numeric slopedness parameter.
#'
#' @return A vector containing a value of \code{level} for every element in
#'   \code{t}.
#' @export
#'
#' @examples
#' unif_traj(0:12)
logistic_traj = function(t, offset=0, a=2)
{
  t = t + offset
  result = rep(0, length(t))
  result[(t %% 12) <= 6] = 10 + 90/(1+exp((3-(t[(t %% 12) <= 6] %% 12)) * a))
  result[(t %% 12) >  6] = 10 + 90/(1+exp(((t[(t %% 12) >  6] %% 12) - 12 + 3) * a))
  return(result)
}

#' @export
logistic_traj_inv = function(t, offset=0, a=2)
{
  return(1/logistic_traj(t, offset, a))
}

#' @export
bottleneck_traj <- function(t)
{
  result = rep(0,length(t))
  result[t<=0.5]<-1
  result[t>0.5 & t<1]<-.1
  result[t>=1]<-1
  return(result)
}

#' @export
bottleneck_traj_inv <- function(t)
{
  return(1/bottleneck_traj(t))
}
