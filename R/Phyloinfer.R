#### Elliptical slice sampler by Murray et~al (2010) ####

# inputs:
#   q_cur: initial state of the parameter
#   l_cur: initial log-likelihood
#   loglik: log-likelihood function of q
#   cholC: Cholesky decomposition (upper triangular matrix) of covariance matrix of Gaussian prior
# outputs:
#   q: new state of the parameter following N(q;0,Cov)*lik
#   l: log-likelihood of new state
#   Ind: proposal acceptance indicator

ESS = function(q_cur, l_cur, loglik, cholC)
{  
  # choose ellipse
  nu = t(cholC)%*%rnorm(length(q_cur))
  
  # log-likelihood threshold
  u = runif(1)
  logy <- l_cur + log(u)
  
  # draw a initial proposal, also defining a bracket
  t = 2*pi*runif(1)
  t_min <- t-2*pi; t_max <- t
  
  while(1)
  {
    q <- q_cur*cos(t) + nu*sin(t)
    l = loglik(q)
    if(l>logy) return(list(q=q,l=l,Ind=1))
    # shrink the bracket and try a new point
    if(t<0) t_min <- t
    else t_max <- t
    t = runif(1,t_min,t_max)
  }
}

#### Metropolis-Adjusted Langevin (MALA) Algorithm ####
# This function generates one sample given previous state.

# inputs:
#   q_cur: initial state of the parameter
#   u_cur, du_cur: initial potential energy and its gradient
#   U:=-log(density(q)), potential function of q, or its gradient
#   eps: step size
# outputs:
#   q: new state of the parameter
#   u, du: new potential energy and its gradient
#   Ind: proposal acceptance indicator

MALA = function (q_cur, u_cur, du_cur, U, eps=.2)
{
  # initialization
  q = q_cur
  D = length(q)
  u = u_cur
  du = du_cur
  
  # sample momentum
  p = rnorm(D)
  
  # calculate current energy
  E_cur = u + sum(p^2)/2
  
  # Make a half step for momentum
  p = p - eps/2 * du
  
  # Make a full step for the position
  q = q + eps * p
  
  du = U(q,T)
  # Make a half step for momentum at the end
  p = p - eps/2 * du
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  u = U(q)
  E_prp = u + sum(p^2)/2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  logAP = -E_prp + E_cur
  
  if( is.finite(logAP)&&(log(runif(1))<min(0,logAP)) ) 
    return (list(q = q, u = u, du = du, Ind = 1))
  else 
    return (list(q = q_cur, u = u_cur, du = du_cur, Ind = 0))
}

#### adaptive Metropolis-Adjusted Langevin (aMALA) Algorithm ####
# This is adaptive block updating GMRF by Knorr-Held and Rue (2002), equivalent to Riemannian MALA by Girolami and Calderhead (2011).
# This function generates one sample given previous state.

# inputs:
#   q_cur: initial state of the parameter
#   u_cur: initial potential energy
#   U:=-log(density(q)), potential function of q, or its gradient
#   Met: Fisher observed(or expected) information matrix of approximating Normal
#   c: parameter to control step size of kappa
#   eps: step size
#   L: number of leapfrogs
# outputs:
#   q: new state of the parameter
#   u: new potential energy
#   Ind: proposal acceptance indicator

aMALA = function (q_cur, u_cur, U, Met, c, eps=1)
{
  # initialization
  q = q_cur
  D = length(q)
  
  # sample kappa
  repeat
  {
    t=runif(1,1/c,c)
    if(runif(1)<(t+1/t)/(c+1/c))
      break
  }
  q[D]=q[D]*t
  
  # prepare pre-conditional matrix and gradient
  Q=Met(q)
  cholQ=chol(Q)
  g=U(q,T)
  
  # sample momentum
  z=rnorm(D-1)
  p=backsolve(cholQ,z)
  
  # log proposal density
  logprp = -t(z)%*%z/2+sum(log(diag(cholQ)))
  
  # update momentum
  #	p=p-eps/2*solve(Q,g[-D])
  p = p-eps/2*(chol2inv(cholQ)%*%g[-D])
  
  # update position
  q[-D] = q[-D]+eps*p
  
  # update pre-conditional matrix and gradient
  Q = Met(c(q[-D],q_cur[D]))
  cholQ=chol(Q)
  g = U(c(q[-D],q_cur[D]),T) # very interesting update!!!
  
  # update momentum
  p = p-eps/2*(chol2inv(cholQ)%*%g[-D])
  
  # log reverse proposal density
  logprp_rev = -t(p)%*%Q%*%p/2+sum(log(diag(cholQ)))
  
  # Evaluate potential energy
  u = U(q)
  
  # Accept or reject the state jointly
  logAP = -u + u_cur - logprp + logprp_rev
  
  if ( is.finite(logAP) && (log(runif(1))<min(0,logAP)) )
    return (list(q = q, u = u, Ind = 1))
  else
    return (list(q = q_cur, u = u_cur, Ind = 0))
}

#### Hamiltonian Monte Carlo ####
# This is standard HMC method.
# This function generates one sample given previous state

# inputs:
#   q_cur: initial state of the parameter
#   u_cur, du_cur: initial potential energy and its gradient
#   U:=-log(density(q)), potential function of q, or its gradient
#   eps: step size
#   L: number of leapfrogs
# outputs:
#   q: new state of the parameter
#   u, du: new potential energy and its gradient
#   Ind: proposal acceptance indicator

HMC = function (q_cur, u_cur, du_cur, U, eps=.2, L=5, rand_leap=TRUE)
{  
  # initialization
  q = q_cur
  D = length(q)
  u = u_cur
  du = du_cur
  
  # sample momentum
  p = rnorm(D)
  
  # calculate current energy
  E_cur = u + sum(p^2)/2
  
  # Make a half step for momentum at the beginning
  p = p - eps/2 * du
  
  if (rand_leap)
    randL = ceiling(runif(1)*L)
  else
    randL = ceiling(L)
  
  # Alternate full steps for position and momentum
  for (l in 1:randL)
  {
    # Make a full step for the position
    q = q + eps * p
    
    du = U(q,T)
    # Make a full step for the momentum, except at end of trajectory
    if (l!=randL)
      p = p - eps * du
  }
  
  # Make a half step for momentum at the end.
  p = p - eps/2 * du
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  u = U(q)
  E_prp = u + sum(p^2)/2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  logAP = -E_prp + E_cur
  
  if( is.finite(logAP) && (log(runif(1))<min(0,logAP)) )
    return (list(q = q, u = u, du = du, Ind = 1))
  else
    return (list(q = q_cur, u = u_cur, du = du_cur, Ind = 0))
}

#### Split Hamiltonian Monte Carlo ####
# This is splitHMC method by (Gaussian) approximation.
# This function generates one sample given previous state.

# inputs:
#   q_cur: initial state of the parameter
#   u_cur, du_cur: initial potential energy and its gradient
#   U:=-log(density(q)), potential function of q, or its gradient
#   rtEV, EVC: square root of eigen-valudes, eigen-vectors of Fisher observed(or expected) information matrix of approximating Normal
#   eps: step size
#   L: number of leapfrogs
# outputs:
#   q: new state of the parameter
#   u, du: new potential energy and its gradient
#   Ind: proposal acceptance indicator

splitHMC = function (q_cur, u_cur, du_cur, U, rtEV, EVC, eps=.1, L=5, rand_leap=TRUE)
{
  # initialization
  q = q_cur
  D = length(q)
  u = u_cur
  du = du_cur
  
  # sample momentum
  p = rnorm(D)
  
  # calculate current energy
  E_cur = u + sum(p^2)/2
  
  
  if (rand_leap)
    randL = ceiling(runif(1)*L)
  else
    randL = ceiling(L)
  
  p = p - eps/2*du
  qT = rtEV*(t(EVC)%*%q[-D])
  pT = t(EVC)%*%p[-D]
  A = t(qT)%*%qT
  # Alternate full steps for position and momentum
  for (l in 1:randL)
  {
    p[D] <- p[D] - eps/2*A/2*exp(q[D])
    q[D] <- q[D] + eps/2*p[D]
    
    # Make a full step for the middle dynamics
    Cpx = complex(mod=1,arg=-rtEV*exp(q[D]/2)*eps)*complex(re=qT*exp(q[D]/2),im=pT)
    qT = Re(Cpx)*exp(-q[D]/2)
    pT = Im(Cpx)
    q[-D] = EVC%*%(qT/rtEV)
    
    # Make a half step for the last half dynamics
    A=t(qT)%*%qT
    
    q[D] <- q[D] + eps/2*p[D]
    p[D] <- p[D] - eps/2*A/2*exp(q[D])
    
    du = U(q,T)
    if(l!=randL)
    {
      pT = pT - eps*(t(EVC)%*%du[-D])
      p[D] = p[D] - eps*du[D]
    }
  }
  p[-D] = EVC%*%pT - eps/2*du[-D]
  p[D] = p[D] - eps/2*du[D]
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  u = U(q)
  E_prp = u + sum(p^2)/2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  logAP = -E_prp + E_cur
  
  if( is.finite(logAP) && (log(runif(1))<min(0,logAP)) )
    return (list(q = q, u = u, du = du, Ind = 1))
  else
    return (list(q = q_cur, u = u_cur, du = du_cur, Ind = 0))
}

#### Helper functions

# Intrinsic precision matrix
Q_matrix <- function(input,s.noise,signal)
{
  n2<-nrow(input)
  diff1<-diff(input)
  diff1[diff1==0]<-s.noise #correction for dividing over 0
  diff<-(1/(signal*diff1))
  Q<-spam(0,n2,n2)  
  if (n2>2)
  {
    Q[cbind(seq(1,n2),seq(1,n2))] <- c(diff[1],diff[1:(n2-2)]+diff[2:(n2-1)],diff[n2-1])+(1/signal)*rep(s.noise,n2)
  }
  else
  {
    Q[cbind(seq(1,n2),seq(1,n2))] <- c(diff[1],diff[n2-1])+(1/signal)*rep(s.noise,n2)
  }
  Q[cbind(seq(1,n2-1),seq(2,n2))] <- -diff[1:(n2-1)]
  Q[cbind(seq(2,n2),seq(1,n2-1))] <- -diff[1:(n2-1)]
  return(Q)
}

# backwards compatibility (deprecate soon)
Q.matrix = function(...)
{
  return(Q_matrix(...))
}

#### Sampling wrappers ####

# This serves as a black box to sample distributions using HMC algorithms provided data and basic settings. #
sampling = function(data, para, alg, setting, init, print=TRUE)
{
  # pass the data and parameters
  lik_init = data$lik_init # f_offset = data$f_offset
  Ngrid = lik_init$ng+1
  alpha = para$alpha
  beta = para$beta
  invC = para$invC
  rtEV = para$rtEV
  EVC = para$EVC
  cholC = para$cholC
  
  # MCMC sampling setting
  stepsz = setting$stepsz
  Nleap  = setting$Nleap
  if (alg=='aMALA')
    szkappa=setting$szkappa
  if (alg=="HMC" | alg == "splitHMC")
    rand_leap = settings$rand_leap
  
  # storage of posterior samples
  NSAMP = setting$NSAMP
  NBURNIN = setting$NBURNIN
  samp = matrix(NA,NSAMP-NBURNIN,Ngrid) # all parameters together
  acpi = 0
  acpt = 0
  
  # initialization
  theta = init$theta
  u = init$u
  du = init$du
  
  # start MCMC run
  start_time = Sys.time()
  cat('Running ', alg ,' sampling...\n')
  for(Iter in 1:NSAMP)
  {  
    if(print&&Iter%%100==0)
    {
      cat(Iter, ' iterations have been finished!\n' )
      cat('Online acceptance rate is ',acpi/100,'\n')
      acpi=0
    }
    
    # sample the whole parameter
    tryCatch({res=switch(alg,
                         HMC=eval(parse(text='HMC'))(theta,u,du,function(theta,grad=FALSE)U(theta,lik_init,invC,alpha,beta,grad),stepsz,Nleap,rand_leap),
                         splitHMC=eval(parse(text='splitHMC'))(theta,u,du,function(theta,grad=FALSE)U_split(theta,lik_init,invC,alpha,beta,grad),rtEV,EVC,stepsz,Nleap,rand_leap),
                         MALA=eval(parse(text='MALA'))(theta,u,du,function(theta,grad=FALSE)U(theta,lik_init,invC,alpha,beta,grad),stepsz),
                         aMALA=eval(parse(text='aMALA'))(theta,u,function(theta,grad=FALSE)U_kappa(theta,lik_init,invC,alpha,beta,grad),function(theta)Met(theta,lik_init,invC),szkappa,stepsz),
                         ESS=eval(parse(text='ESS'))(theta[-Ngrid],u,function(f)coal_loglik(lik_init,f),cholC/sqrt(theta[Ngrid])),
                         stop('The algorithm is not in the list!'));
              theta[1:(Ngrid-(alg=='ESS'))]=res$q;u=res[[2]];if(any(grepl(alg,c('HMC','splitHMC','MALA'))))du=res$du;
              acpi=acpi+res$Ind}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    # Gibbs sample kappa for ESS
    if(alg=='ESS')
      theta[Ngrid]=rgamma(1,alpha+(Ngrid-1)/2,beta+t(theta[-Ngrid])%*%invC%*%theta[-Ngrid]/2)
    
    # save posterior samples after burnin
    if(Iter>NBURNIN)
    {
      samp[Iter-NBURNIN,]<-theta
      acpt<-acpt+res$Ind
    }
    
  }
  stop_time = Sys.time()
  time = stop_time-start_time
  cat('\nTime consumed : ',time)
  acpt = acpt/(NSAMP-NBURNIN)
  cat('\nFinal Acceptance Rate: ',acpt,'\n')
  
  return(list(samp=samp,time=time,acpt=acpt))
}

# Same as 'sampling' above but specially designed for comparing mixing rate of MCMC algorithms #
sampling_mixrate = function(data, para, alg, setting, init, print=TRUE)
{
  # pass the data and parameters
  lik_init = data$lik_init # f_offset = data$f_offset
  Ngrid = lik_init$ng+1
  alpha = para$alpha
  beta = para$beta
  invC = para$invC
  rtEV = para$rtEV
  EVC = para$EVC
  cholC = para$cholC
  
  # MCMC sampling setting
  stepsz = setting$stepsz
  Nleap  = setting$Nleap
  if (alg=='aMALA')
    szkappa = setting$szkappa
  if (alg=="HMC" | alg=="splitHMC")
    rand_leap = setting$rand_leap
  
  # storage of posterior samples
  WallTime = setting$WallTime
  Intvl = setting$Intvl
  SaveLeng = ceiling(WallTime/Intvl)
  acpi = 0
  acpt = 0
  
  # save for comparing mixing rate
  logLiks = rep(NA,SaveLeng) # save all the log-likelihoods
  times = rep(NA,SaveLeng) # save all time intervals
  
  # initialization
  theta = init$theta
  u  = init$u
  du = init$du
  times[1] = 0
  logLiks[1] = coal_loglik(lik_init,theta[-Ngrid])
  
  # start MCMC run
  start_time = Sys.time()
  cat('Running ', alg ,' sampling...\n')
  Iter=1
  counter=1
  while(counter<=SaveLeng&times[counter]<=WallTime)
  {    
    if(print&&Iter%%100==0)
    {
      cat(Iter, ' iterations have been finished!\n' )
      cat('Online acceptance rate is ',acpi/100,'\n')
      acpi=0
    }
    
    # sample the whole parameter
    tryCatch({res=switch(alg,
                         HMC=eval(parse(text='HMC'))(theta,u,du,function(theta,grad=FALSE)U(theta,lik_init,invC,alpha,beta,grad),stepsz,Nleap,rand_leap),
                         splitHMC=eval(parse(text='splitHMC'))(theta,u,du,function(theta,grad=FALSE)U_split(theta,lik_init,invC,alpha,beta,grad),rtEV,EVC,stepsz,Nleap,rand_leap),
                         MALA=eval(parse(text='MALA'))(theta,u,du,function(theta,grad=FALSE)U(theta,lik_init,invC,alpha,beta,grad),stepsz),
                         aMALA=eval(parse(text='aMALA'))(theta,u,function(theta,grad=FALSE)U_kappa(theta,lik_init,invC,alpha,beta,grad),function(theta)Met(theta,lik_init,invC),szkappa,stepsz),
                         ESS=eval(parse(text='ESS'))(theta[-Ngrid],u,function(f)coal_loglik(lik_init,f),cholC/sqrt(theta[Ngrid])),
                         stop('The algorithm is not in the list!'));
              theta[1:(Ngrid-(alg=='ESS'))]=res$q;u=res[[2]];if(any(grepl(alg,c('HMC','splitHMC','MALA'))))du=res$du;
              acpi=acpi+res$Ind}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    # Gibbs sample kappa for ESS
    if(alg=='ESS')
      theta[Ngrid]=rgamma(1,alpha+(Ngrid-1)/2,beta+t(theta[-Ngrid])%*%invC%*%theta[-Ngrid]/2)
    
    # save acceptance rate
    acpt<-acpt+res$Ind
    
    cur_time=as.numeric(Sys.time()-start_time,units='secs')
    if(cur_time-times[counter]>Intvl)
    {
      counter=counter+1
      times[counter]=cur_time
      logLiks[counter]=coal_loglik(lik_init,theta[-Ngrid])
    }
    Iter = Iter+1;
  }
  stop_time = Sys.time()
  time = stop_time-start_time
  cat('\nTime consumed : ',time)
  acpt = acpt/Iter
  cat('\nFinal Acceptance Rate: ',acpt,'\n')
  
  return(list(time=time,acpt=acpt,WallTime=WallTime,Intvl=Intvl,logLiks=logLiks,times=times))
}

# wrapper that encapsulates sampler above with good defaults
mcmc_sampling = function(data, alg, nsamp, nburnin, Ngrid=100, nugget="1,1", prec_alpha = 1e-2, prec_beta = 1e-2,
                         TrjL=NULL, Nleap=NULL, szkappa=NULL, rand_leap=NULL)
{
  # add ability to parse genealogy objects as well as lists
  samp_times = data$samp_times
  n_sampled  = data$n_sampled
  coal_times = data$coal_times
  
  # Jump tuning parameters--should probably have an option to change in the arguments
  if (is.null(TrjL))
    TrjL = switch(alg, HMC=3, splitHMC=3, MALA=0.1, aMALA=0.1)
  if (is.null(Nleap))
    Nleap = switch(alg, HMC=30, splitHMC=15, MALA=1, aMALA=1)
  if (is.null(szkappa) & alg=="aMALA")
    szkappa = 1.2
  
  stepsz = TrjL/Nleap
  
  grid_bds = range(c(coal_times,samp_times))
  #Ngrid = 100
  
  grid = seq(grid_bds[1],grid_bds[2],length.out=Ngrid)
  intl = grid[2]-grid[1]
  midpts = grid[-1]-intl/2
  
  # initialize likelihood calculation
  lik_init = coal_lik_init(samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times, grid=grid)
  
  # calculate intrinsic precision matrix
  invC <- Q.matrix(as.matrix(midpts),0,1)
  
  # fudge to be able to compute the cholC
  if (nugget == "1,1")
    invC[1,1] <- invC[1,1]+.0001 # nugget at (1,1)
  else if (nugget == "diag")
    diag(invC)<-diag(invC)+.0001 # nugget for the whole diagonal
  else if (nugget == "none")
    warning("No nugget may result in a non-full-rank matrix.")
  else
    stop(paste("Unrecognized argument nugget = '", nugget, "', please use '1,1', 'diag', or 'none'.", sep = ""))
  
  eig  = eigen(invC,T)
  rtEV = sqrt(eig$values)
  EVC  = eig$vectors
  
  C = matrix(midpts,Ngrid-1,Ngrid-1)
  C = matrix(pmin(C[col(C)],C[row(C)]),Ngrid-1,Ngrid-1)
  cholC = chol(C)
  
  # initializations
  theta = rep(1,Ngrid)
  u  = U_split(theta,lik_init,invC,prec_alpha,prec_beta)
  du = U_split(theta,lik_init,invC,prec_alpha,prec_beta,TRUE)
  
  # MCMC sampling preparation
  data = list(lik_init=lik_init)
  para = list(alpha=prec_alpha,beta=prec_beta,invC=invC,rtEV=rtEV,EVC=EVC,cholC=cholC)
  setting = data.frame(stepsz=stepsz,Nleap=Nleap,NSAMP=nsamp,NBURNIN=nburnin,rand_leap=rand_leap)
  init = list(theta=theta,u=u,du=du)
  
  # Run MCMC sampler
  res_MCMC = sampling(data,para,alg,setting,init)
  
  # estimates given by MCMC samples
  med = apply(exp(res_MCMC$samp[,-Ngrid]),2,median);
  low = apply(exp(res_MCMC$samp[,-Ngrid]),2,function(x)quantile(x,.025))
  up  = apply(exp(res_MCMC$samp[,-Ngrid]),2,function(x)quantile(x,.975))
  
  # Incorporate estimates into the output from sampling
  res_MCMC$grid = grid
  res_MCMC$med = med
  res_MCMC$low = low
  res_MCMC$up  = up
  res_MCMC$med_fun = stepfun(grid,c(0,med,0))
  res_MCMC$low_fun = stepfun(grid,c(0,low,0))
  res_MCMC$up_fun  = stepfun(grid,c(0,up,0))
  
  return(res_MCMC)
}
