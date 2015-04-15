coalgen_thinning_iso <- function(sample,trajectory,upper=25,...)
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
    if (runif(1)<=trajectory(time,...)/upper)
    {
      out[n-j+1] <- time
      j <- j-1
    }
  }
  return(list(intercoal_times=c(out[1],diff(out)),lineages=seq(n,2,-1), coal_times=out))
}

coalgen_thinning_hetero <- function(sample,trajectory,upper=25,...)
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
    if (runif(1) <= trajectory(E+s,...)/upper)
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
    if (runif(1)<=trajectory(E+s,...)/upper)
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
  
  return(list(intercoal_times=c(out[1],diff(out)),lineages=branches,coal_times=out))   
}
