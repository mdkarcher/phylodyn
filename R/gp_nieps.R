lambda_prior = function(x,lambda_mean,alpha1,alpha2,alpha3)
{
  indicator = rep(1,length(x))
  indicator[x >= lambda_mean] = rep(0,sum(x >= lambda_mean))
  return(alpha2*indicator + alpha3 * exp(-alpha1*x) * (1-indicator))
}

gp_nieps_iso = function(filename, alpha_prior_kappa, beta_prior_kappa, iterations, lambda_mean)
{
  data1 <- utils::read.delim(filename, header=FALSE, sep="\t")
  data <- cbind(cumsum(data1[,1]),data1[,2])
  s <- c(0,data[,1])
  T <- max(s)
  grid.points1 <- seq(0.001,T+.001,length.out=150)
  grid.points2 <- seq(0.001,T+.001,length.out=250)
  grid.points3 <- seq(0.001,T+.001,length.out=500)
  
  params <- matrix(nrow=iterations,ncol=2)
  
  result1 <- matrix(nrow=length(grid.points1),ncol=iterations)
  result2 <- matrix(nrow=length(grid.points2),ncol=iterations)
  result3 <- matrix(nrow=length(grid.points3),ncol=iterations)
  
  loglik <- rep(0,iterations)
  
  alpha2 <- .01/lambda_mean
  alpha3 <- exp(1)*0.99/lambda_mean
  alpha1 <- 1/lambda_mean
  
  
  ## Initial state (consider making arguments)
  s.noise <- .000001
  lambda <- 2
  signal <- lambda/3
  data <- cbind(data,rep(1,nrow(data)),GP.prior(signal,data[,1],s.noise))
  data[,4] <- data[,4]-mean(data[,4])
  l <- length(s)
  m <- rep(0,l-1)
  
  n <- nrow(data)+1
  coal.factor <- rep(0,n-1)
  t <- rep(0,n-1)
  for (k in n:2)
  {
    coal.factor[n-k+1] <- k*(k-1)/2
    t[n-k+1] <- s[n+1-k+1]-s[n+1-k]
  }
  beta.cum <- sum(t*coal.factor)
  
  
  #### MCMC run
  start_time = Sys.time()
  for (Nsim in 1:iterations)
  {
    tdata <- number.thinned(s,T,lambda,signal,data,s.noise,n,t,coal.factor,m)
    data <- tdata$info
    m <- tdata$m
    loglik[Nsim] <- tdata$loglik
    
    m2 <- nrow(data)-n+1
    if (m2>0)
    {
      data <- location.thinned.uniform2(s,T,lambda,signal,n,data,s.noise,t,m)
    }
    
    slice <- slice.sampling(data,signal,s.noise)
    data <- slice$data
    #result2[Nsim,3] <- sum((num.3-data[,4])!=0)
    alpha_post_kappa <- alpha_prior_kappa+nrow(slice$Q)*.5
    gg <- data[,4][slice$order]
    beta_post_kappa <- beta_prior_kappa+.5*signal*t(gg)%*%slice$Q%*%gg
    result1[,Nsim] <- sigmoidal(GP.posterior(as.matrix(data[,1]),data[,4],signal,as.matrix(grid.points1),s.noise)$g)*lambda
    result2[,Nsim] <- sigmoidal(GP.posterior(as.matrix(data[,1]),data[,4],signal,as.matrix(grid.points2),s.noise)$g)*lambda
    result3[,Nsim] <- sigmoidal(GP.posterior(as.matrix(data[,1]),data[,4],signal,as.matrix(grid.points3),s.noise)$g)*lambda
    
    signal <- 1/stats::rgamma(1,alpha_post_kappa,beta_post_kappa)
    params[Nsim,1] <- signal
    
    lambda1 <- stats::runif(1,lambda-.005,lambda+.005)
    if (lambda1<0)
    {
      lambda1 <- -lambda1
    }
    rat1 <- lambda_prior(lambda1, lambda_mean, alpha1, alpha2, alpha3) / lambda.prior(lambda, lambda_mean, alpha1, alpha2, alpha3)
    ac <- stats::runif(1)<rat1*((lambda1/lambda)^(nrow(data))*exp(-(lambda1-lambda)*beta.cum))
    if (ac=="TRUE")
    {
      lambda <- lambda1
    }
    
    params[Nsim,2] <- lambda
  }
  
  stop_time = Sys.time()
  
  return(list(result1, result2, result3, runtime=(stop_time-start_time)))
}