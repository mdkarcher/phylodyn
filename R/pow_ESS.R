log_inorm_prior3 <- function(x, prec, mu=rep(0, length(x)))
{
  x_mu = x - mu
  return(sum(dnorm(x_mu, sd = prec ^ -0.5, log = TRUE)))
}

log_betas_prior3 <- function(betas, betas_prec = 0.01)
{
  return(log_inorm_prior3(x = betas, prec = betas_prec))
}

samp_loglik3 = function(init, logpop, betas_list, covars_list = NULL)
{
  logpop = rep(logpop, init$gridrep)
  betas = betas_list$betas
  beta0 = betas[1]
  beta1 = betas[2]
  
  if (!is.null(covars_list) && !is.null(covars_list$covar_vals))
  {
    covar_vals = as.matrix(covars_list$covar_vals)
    covar_vals = as.matrix(covar_vals[rep(1:init$ng, init$gridrep),])
    covar_betas = betas_list$covar_betas
    covars_bs = covar_vals %*% covar_betas
  }
  else
  {
    covars_bs = 0
  }
  
  if (!is.null(covars_list) && !is.null(covars_list$pow_covar_vals))
  {
    pow_covar_vals = as.matrix(covars_list$pow_covar_vals)
    pow_covar_vals = as.matrix(pow_covar_vals[rep(1:init$ng, init$gridrep),])
    pow_covar_betas = betas_list$pow_covar_betas
    pcovars_bs = (beta1 + pow_covar_vals %*% pow_covar_betas) * logpop
  }
  else
  {
    pcovars_bs = beta1 * logpop
  }
  
  fs_betas = covars_bs + pcovars_bs
  
  #llsampevents = beta1 * init$count * f
  #llsampnoevents = init$D * exp(beta0) * exp(f)^beta1
  #llsamp = init$ns * beta0 + sum(llsampevents[!is.na(init$count)]) - sum(llsampnoevents[!is.na(init$count)])
  
  llsampevents = init$count * fs_betas
  llsampnoevents = init$D * exp(beta0 + fs_betas)
  llsamp = init$ns * beta0 + sum(llsampevents[!is.na(init$count)]) - sum(llsampnoevents[!is.na(init$count)])
  #print(init$ns * beta0)
  #print(sum(llsampevents[!is.na(init$count)]))
  #print(-1 * sum(llsampnoevents[!is.na(init$count)]))
  
  return(llsamp)
}


ESS_betas_ll3 = function(f, lik_init, betas_list, covars_list = NULL) 
{
  return(coal_loglik(init = lik_init, f = f) +
           samp_loglik3(init = lik_init, logpop = f, betas_list = betas_list,
                        covars_list = covars_list))
}

compute_pos_summ3 = function(f, prec, hyperparams, lik_init, setting,
                             betas_list = NULL, covars_list = NULL)
{
  if (setting$samp_alg == "MH" || setting$samp_alg == "fixed")
  {
    loglik = ESS_betas_ll3(f = f, lik_init = lik_init, betas_list = betas_list,
                           covars_list = covars_list)
  }
  else
    loglik = ESS_none_ll(f = f, lik_init = lik_init)
  
  result = data.frame(loglik)
  
  logfieldpri = log_field_prior(f = f, prec = prec, first_elem_prec = hyperparams$first_elem_prec)
  result$logfieldpri = logfieldpri
  
  logpri = logfieldpri
  
  if (!setting$noprec)
  {
    logprecpri = log_kappa_prior(kappa = prec, alpha = hyperparams$alpha, beta = hyperparams$beta)
    result$logprecpri = logprecpri
    
    logpri = logpri + logprecpri
  }
  
  if (setting$samp_alg == "MH")
  {
    logbetapri = log_betas_prior3(betas = c(betas_list$betas,
                                            betas_list$covar_betas,
                                            betas_list$pow_covar_betas),
                                  betas_prec = hyperparams$betas_prec)
    result$logbetapri = logbetapri
    
    logpri = logpri + logbetapri
  }
  
  result$logpri = logpri
  
  logpos = logpri + loglik
  result$logpos = logpos
  
  return(result)
}

MH_betas_rscan3 = function(curr_betas_list, curr_pos_summ, lik_init, f, prec,
                           hyperparams, setting, covars_list = NULL)
{
  curr_betas = c(curr_betas_list$betas, curr_betas_list$covar_betas,
                curr_betas_list$pow_covar_betas)
  new_betas = c(curr_betas_list$betas, curr_betas_list$covar_betas,
                curr_betas_list$pow_covar_betas)
  n_covar = length(curr_betas_list$covar_betas)
  n_pow = length(curr_betas_list$pow_covar_betas)
  n = 2 + n_covar + n_pow
  
  if (is.null(setting$niter))
    niter = n
  
  rscan = sample(x = 1:n, size = niter, replace = TRUE)
  inds = vector(mode="list", length=n)
  
  for (idx in rscan)
  {
    new_betas[idx] = stats::rnorm(n = 1, mean = curr_betas[idx], sd = setting$proposal_sds)
    new_betas_list = split(new_betas,
                           rep(c("betas", "covar_betas",
                                 "pow_covar_betas"),
                               c(2, n_covar, n_pow)))
    
    new_pos_summ = compute_pos_summ3(f = f, prec = prec, hyperparams = hyperparams,
                                     lik_init = lik_init, setting = setting,
                                     betas_list = new_betas_list,
                                     covars_list = covars_list)
    
    if (new_pos_summ$logpos > curr_pos_summ$logpos ||
        log(stats::runif(n = 1)) < new_pos_summ$logpos - curr_pos_summ$logpos)
    {
      curr_betas = new_betas
      curr_pos_summ = new_pos_summ
      inds[[idx]] = c(inds[[idx]], 1)
    }
    else
    {
      new_betas = curr_betas
      inds[[idx]] = c(inds[[idx]], 0)
    }
  }
  
  curr_betas_list = split(curr_betas,
                          rep(c("betas", "covar_betas",
                                "pow_covar_betas"),
                              c(2, n_covar, n_pow)))
  
  result = list(betas_list = curr_betas_list, pos_summ = curr_pos_summ, inds = inds)
  
  return(result)
}


ESS3 = function(q_cur, l_cur, loglik, prec, first_elem_prec, ...)
{  
  # choose ellipse
  # nu = crossprod(cholC, stats::rnorm(length(q_cur)))
  first_elem = rnorm(1, mean = 0, sd = sqrt(1 / first_elem_prec))
  diffs = c(0, rnorm(length(q_cur) - 1, mean = 0, sd = sqrt(1 / prec)))
  nu = first_elem + cumsum(diffs)
  
  # log-likelihood threshold
  u = stats::runif(1)
  logy <- l_cur + log(u)
  
  # draw a initial proposal, also defining a bracket
  t = 2*pi*stats::runif(1)
  t_min <- t-2*pi
  t_max <- t
  
  q <- q_cur*cos(t) + nu*sin(t)
  l <- loglik(q, ...)
  
  while (l < logy)
  {
    # shrink the bracket and try a new point
    if (t < 0)
    {
      t_min <- t
    }
    else
    {
      t_max <- t
    }
    
    t <- stats::runif(1, t_min, t_max)
    q <- q_cur*cos(t) + nu*sin(t)
    
    l <- loglik(q, ...)
  }
  
  return(list(q=q, u=l, Ind=1))
}


sampling_ESS3 = function(data, para, setting, init,
                         verbose=TRUE, printevery=100)
{
  # pass the data and parameters
  lik_init = data$lik_init
  covar_vals = data$covar_vals
  pow_covar_vals = data$pow_covar_vals
  Ngrid = lik_init$ng+1
  
  alpha = para$alpha
  beta = para$beta
  first_elem_prec = para$first_elem_prec # TODO: make mcmc_sampling do this
  beta_prec = para$beta_prec
  
  # storage of posterior samples
  NSAMP = setting$NSAMP
  NBURNIN = setting$NBURNIN
  NSUBSAMP = setting$NSUBSAMP
  recorded_iters = seq.int(from = NBURNIN+1, to = NSAMP, by = NSUBSAMP)
  
  proposal_sds = setting$proposal_sds
  samp_alg = setting$samp_alg
  prec_alg = setting$prec_alg
  printevery = setting$printevery
  
  noprec = FALSE
  setting$noprec = FALSE
  if (prec_alg == "none")
  {
    noprec = TRUE
    setting$noprec = TRUE
  }
  
  # initialization
  f = init$f
  prec = init$prec
  
  fmat = matrix(NA, nrow = length(recorded_iters), ncol = length(f))
  precs = rep(NA, length(recorded_iters))
  
  # storage of log prior, log likelihood, and log posterior
  # logpri = rep(NA, length(recorded_iters))
  # loglik = rep(NA, length(recorded_iters))
  # logpos = rep(NA, length(recorded_iters))
  pos_summ_out = data.frame()
  
  if (samp_alg == "none")
  {
    ll = ESS_none_ll
  }
  else if (samp_alg == "MH" || samp_alg == "fixed")
  {
    betas_list = list(betas = init$betas, covar_betas = init$covar_betas,
                      pow_covar_betas = init$pow_covar_betas)
    ll = ESS_betas_ll3
  }
  
  if (samp_alg == "MH")
  {
    betas_out = matrix(NA, nrow = length(recorded_iters),
                       ncol = 2 + length(betas_list$covar_betas) + length(betas_list$pow_covar_betas))
  }
  
  pos_summ = compute_pos_summ3(f = f, prec = prec, hyperparams = para,
                               lik_init = lik_init, setting = setting,
                               betas_list = betas_list, covars_list = data)
  init_pos = list(f = f, prec = prec, first_elem_prec = first_elem_prec, pos_summ = pos_summ)
  
  # start MCMC run
  start_time = Sys.time()
  cat('Running ESS', samp_alg ,' sampling...\n')
  for(iter in 1:NSAMP)
  {
    if (samp_alg == "none")
    {
      res = ESS3(q_cur = f, l_cur = pos_summ$loglik, loglik = ll, prec = prec,
                 first_elem_prec = first_elem_prec, lik_init = lik_init)
      f = res$q
    }
    else if (samp_alg == "fixed")
    {
      res = ESS3(q_cur = f, l_cur = pos_summ$loglik, loglik = ll, prec = prec,
                 first_elem_prec = first_elem_prec, lik_init = lik_init,
                 betas_list = betas_list, covars_list = data)
      f = res$q
    }
    else if (samp_alg == "MH")
    {
      MH_res = MH_betas_rscan3(curr_betas_list = betas_list,
                               curr_pos_summ = pos_summ,
                               lik_init = lik_init, f = f, prec = prec,
                               hyperparams = para, setting = setting,
                               covars_list = data)
      betas_list = MH_res$betas_list
      pos_summ = MH_res$pos_summ
      
      res = ESS3(q_cur = f, l_cur = pos_summ$loglik, loglik = ll, prec = prec,
                 first_elem_prec = first_elem_prec, lik_init = lik_init,
                 betas_list = betas_list, covars_list = data)
      f = res$q
    }
    else
    {
      stop('The ESS subalgorithm is not in the list!')
    }
    
    #acpi <- acpi+res$Ind
    
    if (prec_alg == "gibbs")
    {
      prec <- stats::rgamma(1, alpha + (Ngrid-1)/2, beta + sum(diff(f)^2)/2)
    }
    else if (prec_alg != "none")
    {
      stop("Precision operator not recognized.")
    }
    
    pos_summ = compute_pos_summ3(f = f, prec = prec, hyperparams = para,
                                 lik_init = lik_init, setting = setting,
                                 betas_list = betas_list, covars_list = data)
    
    # save posterior samples after burnin
    output_index = match(x = iter, table = recorded_iters)
    if (!is.na(output_index))
    {
      fmat[output_index, ] <- f
      precs[output_index] <- prec
      
      if (samp_alg == "MH")
      {
        betas_out[output_index, ] <- c(betas_list$betas, betas_list$covar_betas,
                                       betas_list$pow_covar_betas)
      }
      #acpt[output_index] <- res$Ind
      
      # logpri[output_index] <- pos_summ$logpri
      # loglik[output_index] <- pos_summ$loglik
      # logpos[output_index] <- pos_summ$logpos
      pos_summ_out = rbind(pos_summ_out, pos_summ)
    }
    
    if (verbose && iter %% printevery == 0)
    {
      cat(iter, ' iterations have been finished!\n' )
      # cat('Online acceptance rate is ', 1,'\n')
      #acpi=0
    }
  }
  stop_time <- Sys.time()
  time <- stop_time-start_time
  cat('\nTime consumed : ', time, units(time))
  #acpt <- acpt/(NSAMP-NBURNIN)
  # cat('\nFinal Acceptance Rate: ', sum(acpt) / length(recorded_iters),'\n')
  
  # pos_summ = data.frame(acpt=acpt, logpri = logpri, loglik = loglik, logpos = logpos)
  
  colnames(fmat) = paste("f", 1:(Ngrid-1), sep = "")
  if (samp_alg == "MH")
  {
    samp = data.frame(fmat, precs, betas_out)
    colnames(samp) = c(paste("f", 1:length(f), sep = ""), "prec", paste("beta", 0:(dim(betas_out)[2]-1)))
    
    params_out = data.frame(precs, betas_out)
    colnames(params_out) = c("prec", paste("beta", 0:(dim(betas_out)[2]-1)))
  }
  else
  {
    samp = data.frame(fmat, precs)
    colnames(samp) = c(paste("f", 1:(Ngrid-1), sep = ""), "prec")
    
    params_out = data.frame(precs)
    colnames(params_out) = "prec"
  }
  
  return(list(latent = fmat, params = params_out, precs = precs, samp=samp,
              alg="ESS3", time=time, pos_summ = pos_summ_out,
              samp_alg = samp_alg, prec_alg = prec_alg, init_pos = init_pos))
}


mcmc_ESS = function(dataset, nsamp, nburnin=0, nsubsamp=1, ngrid=100, printevery=100,
                    f_init = rep(1, ngrid-1), prec = 1, betas=rep(0, 2),
                    covariates=NULL, covar_betas=NULL,
                    power_covariates=NULL, pow_covar_betas=NULL,
                    samp_alg = "none", prec_alg = "gibbs", beta_prec = 0.01,
                    first_elem_prec = 0.01, prec_alpha = 0.01, prec_beta = 0.01)
{
  if (class(dataset) == "phylo")
  {
    phy <- summarize_phylo(dataset)
  }
  else if (all(c("coal_times", "samp_times", "n_sampled") %in% names(dataset)))
  {
    phy <- with(dataset, list(samp_times = samp_times, coal_times = coal_times,
                              n_sampled = n_sampled))
  }
  
  samp_times = phy$samp_times
  n_sampled  = phy$n_sampled
  coal_times = phy$coal_times
  
  grid_bds = range(c(coal_times,samp_times))
  grid = seq(grid_bds[1],grid_bds[2],length.out=ngrid)
  intl = grid[2]-grid[1]
  midpts = grid[-1]-intl/2
  
  covar_vals = NULL
  pow_covar_vals = NULL
  if (!is.null(covariates))
  {
    if (is.null(covar_betas))
      covar_betas = rep(0, length(covariates))
    
    for (fcn in covariates)
    {
      covar_vals = cbind(covar_vals, log(fcn(midpts)), deparse.level = 0)
    }
  }
  if (!is.null(power_covariates))
  {
    if (is.null(pow_covar_betas))
      pow_covar_betas = rep(0, length(power_covariates))
    
    for (fcn in power_covariates)
    {
      pow_covar_vals = cbind(pow_covar_vals, fcn(midpts), deparse.level = 0)
    }
  }
  
  # initialize likelihood calculation
  lik_init = coal_lik_init(samp_times = samp_times, n_sampled = n_sampled,
                           coal_times = coal_times, grid = grid)
  
  # MCMC sampling preparation
  dataset = list(lik_init = lik_init, covar_vals = covar_vals, pow_covar_vals = pow_covar_vals)
  para = list(alpha = prec_alpha, beta = prec_beta, beta_prec = beta_prec,
              first_elem_prec = first_elem_prec)
  setting = list(NSAMP = nsamp, NBURNIN = nburnin, NSUBSAMP = nsubsamp,
                 proposal_sds = 0.3, samp_alg = samp_alg, prec_alg = prec_alg,
                 printevery = printevery)
  init = list(f = f_init, prec = prec, betas = betas,
              covar_betas = covar_betas, pow_covar_betas = pow_covar_betas)
  
  # Run MCMC sampler
  res_MCMC = sampling_ESS3(data = dataset, para = para, setting = setting, init = init)
  
  res_MCMC$alg = "ESS3"
  res_MCMC$samp_alg = samp_alg
  res_MCMC$prec_alg = prec_alg
  res_MCMC$Ngrid = ngrid
  
  #cleaned_res = burnin_subsample(res = res_MCMC, burnin = 0)
  
  logfmat = res_MCMC$latent
  params = res_MCMC$params
  estimates = calculate_estimates(logfmat = logfmat, params = params, grid = grid)
  
  #res_MCMC$cleaned_res = cleaned_res
  res_MCMC$estimates = estimates
  
  res_MCMC$med = estimates$fmed
  res_MCMC$low = estimates$flow
  res_MCMC$hi = estimates$fhi
  
  res_MCMC$med_fun = estimates$fmed_fun
  res_MCMC$low_fun = estimates$flow_fun
  res_MCMC$hi_fun = estimates$fhi_fun
  
  res_MCMC$grid = grid
  res_MCMC$x = midpts
  res_MCMC$samp_times = samp_times
  res_MCMC$n_sampled = n_sampled
  res_MCMC$coal_times = coal_times
  
  return(res_MCMC)
}
