gen_INLA_args <- function(samp_times, n_sampled, coal_times)
{
  if (sum(n_sampled) != length(coal_times) + 1)
    stop("Number sampled not equal to number of coalescent events + 1.")
  
  if (length(intersect(coal_times, samp_times)) > 0)
    warning("Coincident sampling event and coalescent event: results may be unpredictable.")
  
  l <- length(samp_times)
  m <- length(coal_times)
  sorting <- sort(c(samp_times, coal_times), index.return=TRUE)
  
  lineage_change <- c(n_sampled, rep(-1, m))[sorting$ix]
  lineages <- head(cumsum(lineage_change), -1) # remove entry for the post-final-coalescent-event open interval
  coal_factor <- lineages*(lineages-1)/2
  
  event <- c(rep(0, l), rep(1, m))[sorting$ix]
  
  return(list(coal_factor=coal_factor, s=sorting$x, event=event, lineages=lineages))
}

gen_summary = function(coal_times, samp_times, n_sampled)
{
  args = gen_INLA_args(coal_times, samp_times, n_sampled)
  n = length(args$s)
  return(data.frame(cbind(lineages=args$lineages, start_time=args$s[1:(n-1)], stop_time=args$s[2:n], end_event=args$event[2:n], change=diff(c(args$indicator,1)))))
}

#' Bayesian nonparametric phylodynamic reconstruction.
#' 
#' @param data \code{phylo} object or list containing vectors of coalescent
#'   times \code{coal_times}, sampling times \code{samp_times}, and number sampled
#'   per sampling time \code{n_sampled}.
#' @param lengthout numeric specifying number of grid points - 1.
#' @param pref logical. Should the preferential sampling model be used?
#' @param prec_alpha,prec_beta numerics specifying gamma prior for precision 
#'   \eqn{\tau}.
#' @param beta1_prec numeric specifying precision for normal prior on 
#'   \eqn{\beta_1}.
#' @param log_zero numeric specifying approximate log of zero.
#'   
#' @return Phylodynamic reconstruction of (log) effective population size at the
#'   grid points.
#' @export
#'   
#' @examples
#' data("NY_flu")
#' res = BNPR(NY_flu)
#' plot_BNPR(res)
BNPR <- function(data, lengthout = 100, pref=FALSE, prec_alpha=0.01,
                 prec_beta=0.01, beta1_prec = 0.001, simplify = TRUE)
{
  if (class(data) == "phylo")
  {
    phy <- summarize_phylo(data)
    samp_times <- phy$samp_times
    n_sampled  <- phy$n_sampled
    coal_times <- phy$coal_times
  }
  else if (all(c("coal_times", "samp_times", "n_sampled") %in% names(data)))
  {
    samp_times <- data$samp_times
    n_sampled  <- data$n_sampled
    coal_times <- data$coal_times
  }
  
  if (!pref)
  {
    result <- infer_coal(samp_times = samp_times, coal_times = coal_times,
                         n_sampled = n_sampled, lengthout = lengthout,
                         prec_alpha = prec_alpha, prec_beta = prec_beta,
                         simplify = simplify)
    result$effpop    <- exp(-result$result$summary.random$time$`0.5quant`)
    result$effpop975 <- exp(-result$result$summary.random$time$`0.025quant`)
    result$effpop025 <- exp(-result$result$summary.random$time$`0.975quant`)
  }
  else
  {
    result = infer_coal_samp(samp_times = samp_times, coal_times = coal_times,
                             n_sampled = n_sampled, lengthout = lengthout,
                             prec_alpha = prec_alpha, prec_beta = prec_beta,
                             beta1_prec = beta1_prec, simplify = simplify)
    result$effpop    <- exp(-result$result$summary.random$ii$`0.5quant`)
    result$effpop975 <- exp(-result$result$summary.random$ii$`0.025quant`)
    result$effpop025 <- exp(-result$result$summary.random$ii$`0.975quant`)
  }
  
  result$samp_times <- samp_times
  result$n_sampled  <- n_sampled
  result$coal_times <- coal_times
  
  return(result)
}

#' @describeIn BNPR Uses preferential sampling model.
#' @export
BNPR_PS <- function(data, lengthout = 100, prec_alpha=0.01, prec_beta=0.01,
                    beta1_prec = 0.001, simplify = TRUE)
{
  return(BNPR(data = data, lengthout = lengthout, pref = TRUE,
              prec_alpha = prec_alpha, prec_beta = prec_beta,
              beta1_prec = beta1_prec, simplify = simplify))
}

firstloop <- function(coal_factor, s, event, lengthout)
{
  grid <- seq(min(s),max(s),length.out=lengthout+1)
  u <- diff(grid)
  field <- grid[-1]-u/2
  sgrid <- grid
  event_new <- 0
  time <- 0
  where <- 1
  E.factor <- 0
  for (j in 1:lengthout)
  {
    count <- sum(s>sgrid[j] & s<=sgrid[j+1]) 
    if (count>1)
    {
      points <- s[s>sgrid[j] & s<=sgrid[j+1]]
      u <- diff(c(sgrid[j],points))
      event_new <- c(event_new,event[(where):(where+count-1)])
      time <- c(time,rep(field[j],count))
      E.factor <- c(E.factor,coal_factor[where:(where+count-1)]*u)
      where <- where+count
      if (max(points)<sgrid[j+1])
      {
        event_new <- c(event_new,0)
        time <- c(time,field[j])
        E.factor <- c(E.factor,coal_factor[where]*(sgrid[j+1]-max(points)))
      }
    }
    if (count==1)
    {
      event_new <- c(event_new,event[where])
      points <- s[s>sgrid[j] & s<=sgrid[j+1]]
      if (points==sgrid[j+1])
      {
        E.factor <- c(E.factor,coal_factor[where]*(sgrid[j+1]-sgrid[j]))
        time <- c(time,field[j])
        where <- where+1
      }
      else
      {
        event_new <- c(event_new,0)
        E.factor <- c(E.factor,coal_factor[where]*(points-sgrid[j]))
        E.factor <- c(E.factor,coal_factor[where+1]*(sgrid[j+1]-points))
        time <- c(time,rep(field[j],2))
        where <- where+1
      }
    }
    if (count==0)
    {
      event_new <- c(event_new,0)
      E.factor <- c(E.factor,coal_factor[where]*(sgrid[j+1]-sgrid[j]))
      time <- c(time,field[j])
    }
  }
  return(list(time=time, event_new=event_new, E_factor=E.factor, grid=grid, field=field))
}

secondloop <- function(time2, event_new2, E_factor2, lengthout, field)
{
  for (j in 1:lengthout)
  {
    count <- sum(time2==field[j])
    if (count>1)
    {
      indic <- seq(1:length(event_new2))[time2==field[j]]
      if (sum(event_new2[indic])==0)
      {
        event_new2 <- event_new2[-indic[-1]]
        time2 <- time2[-indic[-1]]
        temp <- sum(E_factor2[indic])
        E_factor2[indic[1]] <- temp
        E_factor2 <- E_factor2[-indic[-1]]
      }
      #else {}
    }
  }
  
  return(list(time2=time2, event_new2=event_new2, E_factor2=E_factor2))
}

calculate_moller_hetero <- function(coal.factor, s, event, lengthout,
                                    prec_alpha = 0.01, prec_beta = 0.01,
                                    log_zero = -100, alpha = NULL, beta = NULL)
{
  if (prec_alpha == 0.01 & prec_beta==0.01 & !is.null(alpha) & !is.null(beta))
  {
    prec_alpha <- alpha
    prec_beta  <- beta
  }
  
  fl <- firstloop(coal_factor = coal.factor, s = s, event = event,
                  lengthout = lengthout)
  
  sl <- secondloop(time2 = fl$time, event_new2 = fl$event_new,
                   E_factor2 = fl$E_factor, lengthout = lengthout,
                   field = fl$field)
  
  E_log = log(sl$E_factor2)
  E_log[sl$E_factor2 == 0] = log_zero
  
  data <- list(y = sl$event_new2[-1], time = sl$time2[-1], E = E_log[-1])
  formula <- y ~ -1 + f(time, model="rw1", hyper = list(prec = list(param = c(prec_alpha, prec_beta))), constr = FALSE)
  mod4 <- INLA::inla(formula, family = "poisson", data = data, offset = E, control.predictor = list(compute=TRUE))
  
  return(list(result=mod4,grid=fl$grid,data=data,E=E_log, x = fl$field))
}

coal_stats <- function(grid, samp_times, coal_times, n_sampled = NULL,
                       log_zero = -100)
{
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2
  
  if (is.null(n_sampled))
    n_sampled <- rep(1, length(samp_times))
  args <- gen_INLA_args(samp_times = samp_times, n_sampled = n_sampled,
                        coal_times = coal_times)
  
  coal_factor <- args$coal_factor
  s <- args$s
  event <- args$event
  
  grid_trimmed <- setdiff(x = grid, y = s)
  sorting <- sort(c(grid_trimmed, s), index.return=TRUE)
  sgrid <- sorting$x
  ordering <- sorting$ix
  
  time_index <- cut(x = sgrid[-1], breaks = grid, labels = FALSE)
  time <- field[time_index]
  
  event_out <- c(rep(0, length(grid_trimmed)), event)[ordering]
  
  Cfun <- stepfun(x = s, y = c(0, coal_factor, 0), right = TRUE)
  Cvec <- Cfun(sgrid[-1])
  E <- diff(sgrid)*Cvec
  
  E_log = log(E)
  E_log[E == 0] = log_zero
  
  return(data.frame(time = time, event = event_out[-1], E = E, E_log = E_log))
}

condense_stats <- function(time, event, E, log_zero = -100)
{
  result <- aggregate(event ~ time, FUN = sum)
  result$E <- aggregate(E ~ time, FUN = sum)$E
  
  E_log = log(result$E)
  E_log[result$E == 0] = log_zero
  result$E_log <- E_log
  
  return(result)
}

#' @export
infer_coal <- function(samp_times, coal_times, n_sampled = NULL, lengthout = 100,
                       prec_alpha = 0.01, prec_beta = 0.01, simplify = FALSE)
{
  if (min(coal_times) < min(samp_times))
    stop("First coalescent time occurs before first sampling time")
  
  if (max(samp_times) > max(coal_times))
    stop("Last sampling time occurs after last coalescent time")
  
  grid <- seq(min(samp_times), max(coal_times), length.out = lengthout+1)
  
  if (is.null(n_sampled))
    n_sampled <- rep(1, length(samp_times))
  
  coal_data <- coal_stats(grid = grid, samp_times = samp_times, n_sampled = n_sampled,
                          coal_times = coal_times)
  
  if (simplify)
    coal_data <- with(coal_data,
                      condense_stats(time = time, event = event, E=E))
  
  data <- with(coal_data, data.frame(y = event, time = time, E_log = E_log))
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  formula <- y ~ -1 + f(time, model="rw1", hyper = hyper, constr = FALSE)
  
  mod <- INLA::inla(formula, family = "poisson", data = data, offset = data$E_log,
                    control.predictor = list(compute=TRUE))
  
  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}

samp_stats <- function(grid, samp_times, n_sampled = NULL, trim_end = FALSE)
{
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2
  E=diff(grid)
  
  buckets <- cut(x = samp_times, breaks = grid,
                 include.lowest = TRUE)
  
  if (is.null(n_sampled))
    count <- as.vector(table(buckets))
  else
  {
    tab <- aggregate(n_sampled ~ buckets, FUN = sum, labels = FALSE)
    count <- rep(0, lengthout)
    count[as.numeric(tab$buckets)] <- tab$n_sampled
  }
  
  count[head(grid, -1) >= max(samp_times)] <- NA
  result <- data.frame(time = field, count = count, E = E, E_log = log(E))
  
  if (trim_end)
    result <- result[complete.cases(result),]
  
  return(result)
}

#' @export
infer_samp <- function(samp_times, n_sampled = NULL, lengthout = 100,
                       prec_alpha = 0.01, prec_beta = 0.01)
{
  grid <- seq(min(samp_times),max(samp_times),length.out=lengthout+1)
  
  samp_data <- samp_stats(grid = grid, samp_times = samp_times,
                          n_sampled = n_sampled)
  
  data <- with(samp_data, data.frame(y = count, time = time, E_log = E_log))
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  formula_sampling <- y ~ 1 + f(time, model="rw1", hyper = hyper, constr=FALSE)
  
  mod <- INLA::inla(formula_sampling, family="poisson", data=data,
                             offset=data$E_log, control.predictor=list(compute=TRUE))
  
  return(list(result = mod, data = data, grid = grid, x = samp_data$time))
}

calculate_pref <- function(coal.factor,s,event,lengthout,prec_alpha=0.01,prec_beta=0.01)
{
  grid <- seq(0,max(s),length.out=lengthout+1)
  u <- diff(grid)
  field <- grid[-1]-u/2
#   sgrid <- grid
#   event_new <- 0
#   time <- 0
#   where <- 1
#   E.factor <- 0
#   for (j in 1:lengthout)
#   {
#     count <- sum(s>sgrid[j] & s<=sgrid[j+1]) 
#     if (count>1)
#     {
#       points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#       u <- diff(c(sgrid[j],points))
#       event_new <- c(event_new,event[(where):(where+count-1)])
#       time <- c(time,rep(field[j],count))
#       E.factor <- c(E.factor,coal.factor[where:(where+count-1)]*u)
#       where <- where+count
#       if (max(points)<sgrid[j+1])
#       {
#         event_new <- c(event_new,0)
#         time <- c(time,field[j])
#         E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-max(points)))
#       }
#     }
#     if (count==1)
#     {
#       event_new <- c(event_new,event[where])
#       points <- s[s>sgrid[j] & s<=sgrid[j+1]]
#       if (points==sgrid[j+1])
#       {
#         E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
#         time <- c(time,field[j])
#         where <- where+1
#       }
#       else
#       {
#         event_new <- c(event_new,0)
#         E.factor <- c(E.factor,coal.factor[where]*(points-sgrid[j]))
#         E.factor <- c(E.factor,coal.factor[where+1]*(sgrid[j+1]-points))
#         time <- c(time,rep(field[j],2))
#         where <- where+1
#       }
#     }
#     if (count==0)
#     {
#       event_new <- c(event_new,0)
#       E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
#       time <- c(time,field[j])
#     }
#     
#   }
#   time2 <- time
#   event_new2 <- event_new
#   E.factor2 <- E.factor
#   
#   for (j in 1:lengthout)
#   {
#     count <- sum(time2==field[j])
#     if (count>1)
#     {
#       indic <- seq(1:length(event_new2))[time2==field[j]]
#       if (sum(event_new2[indic])==0)
#       {
#         event_new2 <- event_new2[-indic[-1]]
#         time2 <- time2[-indic[-1]]
#         temp <- sum(E.factor2[indic])
#         E.factor2[indic[1]] <- temp
#         E.factor2 <- E.factor2[-indic[-1]]
#       }
#       #else {}
#     }
#   }
  
  #   E.factor2[1:34] <- rep(1,34)
  #data <- list(y=event_new2[-1],event=event_new2[-1],time=time2[-1],E=log(E.factor2[-1]))
  #formula <- y~-1+f(time,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
  #mod4 <- inla(formula,family="poisson",data=data,offset=E,control.predictor=list(compute=TRUE),...)
  
#   n1 <- length(event_new2[-1])
#   n2 <- length(field)
#   Y  <- matrix(NA,n1+n2,2)
#   
  #dd <- heterochronous.gp.stat(influenza.tree)
  #s.time <- dd$sample.times
  #n.sample <- dd$sampled.lineages
  newcount <- rep(0,length(grid)-1)
  
  # MK: added to replace reading from influenza.tree above
  samps = s[event==0]
  
  for (j in 1:length(newcount))
  {
    # MK: altered to samps version
    newcount[j] <- sum(samps > grid[j] & samps <= grid[j+1])
  }
  newcount[1] <- newcount[1]+1
  #newcount[(sum(grid<max(samps))+1):lengthout] <- NA
  newcount[head(grid, -1) >= max(samps)] <- NA
  
  data.sampling2<-data.frame(y=newcount,time=field,E=log(diff(grid)))
  formula.sampling2=y~1+f(time,model="rw1",hyper = list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
  
  # MK: Added during functionizing
  E = diff(s) * coal.factor
  E[1]=1
  
  mod.sampling2 <- INLA::inla(formula.sampling2,family="poisson",data=data.sampling2,offset=E,control.predictor=list(compute=TRUE))
  
  return(list(result=mod.sampling2,grid=grid, data = data.sampling2))
}

joint_stats <- function(coal_data, samp_data)
{
  n1 <- length(coal_data$time)
  n2 <- length(samp_data$time)
  beta0 <- c(rep(0, n1), rep(1, n2))
  E_log <- c(coal_data$E_log, samp_data$E_log)
  Y <- matrix(c(coal_data$event, rep(NA, n2), rep(NA, n1), samp_data$count),
              nrow = n1 + n2, byrow = FALSE)
  w <- c(rep(1, n1), rep(-1, n2))
  ii <- c(coal_data$time, rep(NA, n2))
  jj <- c(rep(NA, n1), samp_data$time)
  
  return(list(Y = Y, beta0 = beta0, ii = ii, jj = jj, w = w, E_log = E_log))
}

#' @export
infer_coal_samp <- function(samp_times, coal_times, n_sampled=NULL, lengthout=100,
                            prec_alpha=0.01, prec_beta=0.01, beta1_prec=0.001,
                            simplify = FALSE, events_only = FALSE)
{
  if (min(coal_times) < min(samp_times))
    stop("First coalescent time occurs before first sampling time")
  
  if (max(samp_times) > max(coal_times))
    stop("Last sampling time occurs after last coalescent time")
  
  grid <- seq(min(samp_times), max(coal_times), length.out = lengthout+1)
  
  if (is.null(n_sampled))
    n_sampled <- rep(1, length(samp_times))
  
  coal_data <- coal_stats(grid = grid, samp_times = samp_times, n_sampled = n_sampled,
                          coal_times = coal_times)
  
  if (simplify)
    coal_data <- with(coal_data, condense_stats(time=time, event=event, E=E))
  
  if (events_only)
    samp_data <- samp_stats(grid = grid, samp_times = samp_times)
  else
    samp_data <- samp_stats(grid=grid, samp_times=samp_times, n_sampled=n_sampled)
  
  joint_data <- joint_stats(coal_data = coal_data, samp_data = samp_data)
  
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  formula <- Y ~ -1 + beta0 +
                 f(ii, model="rw1", hyper = hyper, constr = FALSE) +
                 f(jj, w, copy="ii", fixed=FALSE, param=c(0, beta1_prec))
  mod <- INLA::inla(formula, family = c("poisson", "poisson"),
                    data = joint_data, offset = joint_data$E_log,
                    control.predictor = list(compute=TRUE))
  
  return(list(result = mod, data = joint_data, grid = grid, x = coal_data$time))
}

calculate_moller_hetero_pref <- function(coal.factor,s,event,lengthout,prec_alpha=0.01,prec_beta=0.01,beta1_prec = 0.001,log_zero=-100,alpha=NULL,beta=NULL)
{
  if (prec_alpha == 0.01 & prec_beta==0.01 & !is.null(alpha) & !is.null(beta))
  {
    prec_alpha = alpha
    prec_beta  = beta
  }
  
  grid <- seq(0,max(s),length.out=lengthout+1)
  u <- diff(grid)
  field <- grid[-1]-u/2
  sgrid <- grid
  event_new <- 0
  time <- 0
  where <- 1
  E.factor <- 0
  for (j in 1:lengthout)
  {
    count <- sum(s>sgrid[j] & s<=sgrid[j+1]) 
    if (count>1)
    {
      points <- s[s>sgrid[j] & s<=sgrid[j+1]]
      u <- diff(c(sgrid[j],points))
      event_new <- c(event_new,event[(where):(where+count-1)])
      time <- c(time,rep(field[j],count))
      E.factor <- c(E.factor,coal.factor[where:(where+count-1)]*u)
      where <- where+count
      if (max(points)<sgrid[j+1])
      {
        event_new <- c(event_new,0)
        time <- c(time,field[j])
        E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-max(points)))
      }
    }
    if (count==1)
    {
      event_new <- c(event_new,event[where])
      points <- s[s>sgrid[j] & s<=sgrid[j+1]]
      if (points==sgrid[j+1])
      {
        E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
        time <- c(time,field[j])
        where <- where+1
      }
      else
      {
        event_new <- c(event_new,0)
        E.factor <- c(E.factor,coal.factor[where]*(points-sgrid[j]))
        E.factor <- c(E.factor,coal.factor[where+1]*(sgrid[j+1]-points))
        time <- c(time,rep(field[j],2))
        where <- where+1
      }
    }
    if (count==0)
    {
      event_new <- c(event_new,0)
      E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
      time <- c(time,field[j])
    }
    
  }
  time2 <- time
  event_new2 <- event_new
  E.factor2 <- E.factor
  
  for (j in 1:lengthout)
  {
    count <- sum(time2==field[j])
    if (count>1)
    {
      indic <- seq(1:length(event_new2))[time2==field[j]]
      if (sum(event_new2[indic])==0)
      {
        event_new2 <- event_new2[-indic[-1]]
        time2 <- time2[-indic[-1]]
        temp <- sum(E.factor2[indic])
        E.factor2[indic[1]] <- temp
        E.factor2 <- E.factor2[-indic[-1]]
      }
      #else {}
    }
  }
  
  #   E.factor2[1:34] <- rep(1,34)
  #data <- list(y=event_new2[-1],event=event_new2[-1],time=time2[-1],E=log(E.factor2[-1]))
  #formula <- y~-1+f(time,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
  #mod4 <- inla(formula,family="poisson",data=data,offset=E,control.predictor=list(compute=TRUE),...)
  
  n1 <- length(event_new2[-1])
  n2 <- length(field)
  Y  <- matrix(NA,n1+n2,2)
  
  #dd <- heterochronous.gp.stat(influenza.tree)
  #s.time <- dd$sample.times
  #n.sample <- dd$sampled.lineages
  newcount <- rep(0,length(grid)-1)
  
  # MK: added to replace reading from influenza.tree above
  samps = s[-1][event==0]
  
  for (j in 1:length(newcount))
  {
    # MK: altered to samps version
    newcount[j] <- sum(samps > grid[j] & samps <= grid[j+1])
  }
  newcount[1] <- newcount[1]+1
  newcount[(sum(grid<max(samps))+1):lengthout] <- NA
  
  #data.sampling2<-data.frame(y=newcount,time=field,E=log(diff(grid)))
  #formula.sampling2=y~1+f(time,model="rw1",hyper = list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
  
  # MK: Added during functionizing
  E = diff(s) * coal.factor
  E[1]=1
  
  #print("Got here 1")
  
  beta0<-c(rep(0,n1),rep(1,n2))
  newE<-rep(NA,n1+n2)
  
  #MK: added to resolve INLA failure
  E.factor2.log = log(E.factor2)
  E.factor2.log[E.factor2 == 0] = log_zero
  newE[1:n1]<-E.factor2.log[-1]
  #newE[1:n1]<-log(E.factor2[-1])
  
  newE[(n1+1):(n1+n2)]<-log(diff(grid))
  megafield<-c(time2[-1],field)
  
  #print("Got here 2")
  
  Y[1:n1,1]<-event_new2[-1]
  #print("Got here 3")
  Y[(n1+1):(n2+n1),2]<-newcount
  #print("Got here 4")
  r<-c(rep(1,n1),rep(2,n2))
  w<-c(rep(1,n1),rep(-1,n2))
  # data.pref<-data.frame(Y=Y,beta0=beta0,r=r,megafield=megafield,E=newE,w=w)
  # formula.pref.rep<-Y~-1+beta0+f(megafield,w,model="rw1",replicate=r,hyper=list(prec = list(param = c(.001, .001))),constr=FALSE)
  ii<-c(time2[-1],rep(NA,n2))
  jj<-c(rep(NA,n1),field)
  
  #print("Got here 5")
  
  formula.pref<-Y~-1+beta0+f(ii,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)+
    f(jj,w,copy="ii",fixed=FALSE,param=c(0,beta1_prec))
  #print("Got here 6")
  
  # MK: Changed data.frame to list
  data.pref<-list(Y=Y,beta0=beta0,ii=ii,jj=jj,E=newE,w=w)
  
  #print("Got here 7")
  #print(data.pref)
  mod.pref<-INLA::inla(formula.pref,family=c("poisson","poisson"),offset=E,data=data.pref,control.predictor=list(compute=TRUE))
  
  #print("Got here 8")
  
  return(list(result=mod.pref,data = data.pref, grid=grid, x = field))
}

calculate.moller.hetero.pref = function(...)
{
  return(calculate_moller_hetero_pref(...))
}

gen_INLA_args_old = function(coal_times, s_times, n_sampled)
{
  n         = length(coal_times) + 1
  data      = matrix(0, nrow=n-1, ncol=2)
  data[,1]  = coal_times
  s_times   = c(s_times,max(data[,1])+1)
  data[1,2] = sum(n_sampled[s_times <= data[1,1]])
  tt = length(s_times[s_times <= data[1,1]]) + 1
  
  for (j in 2:nrow(data))
  {
    if (data[j,1] < s_times[tt])
    {
      data[j,2] = data[j-1,2]-1
    }
    else
    {
      data[j,2] = data[j-1,2] - 1 + sum(n_sampled[s_times > data[j-1,1] & s_times <= data[j,1]])
      tt = length(s_times[s_times <= data[j,1]]) + 1
    }  
  }
  
  s = unique(sort(c(data[,1], s_times[1:length(s_times)-1])))
  event1 = sort(c(data[,1], s_times[1:length(s_times)-1]), index.return=TRUE)$ix
  n = nrow(data)+1
  l = length(s)
  event = rep(0, l)
  event[event1<n] = 1
  
  y = diff(s)
  
  coal.factor = rep(0,l-1)
  #indicator = rep(0,l-1) # redundant
  
  t = rep(0,l-1)
  indicator = cumsum(n_sampled[s_times<data[1,1]])
  indicator = c(indicator, indicator[length(indicator)]-1)
  ini = length(indicator) + 1
  for (k in ini:(l-1))
  {
    j = data[data[,1]<s[k+1] & data[,1]>=s[k],2]
    if (length(j) == 0)
    {
      indicator[k] = indicator[k-1] + sum(n_sampled[s_times < s[k+1] & s_times >= s[k]])
    }
    if (length(j) > 0)
    {
      indicator[k] = j-1+sum(n_sampled[s_times < s[k+1] & s_times >= s[k]])
    }
  }
  coal_factor = indicator*(indicator-1)/2
  
  
  return(list(coal_factor=coal_factor, s=s, event=event, lineages=indicator))
}

calculate_moller_hetero_old <- function(coal.factor,s,event,lengthout,prec_alpha=0.01,prec_beta=0.01,log_zero=-100,alpha=NULL,beta=NULL)
{
  if (prec_alpha == 0.01 & prec_beta==0.01 & !is.null(alpha) & !is.null(beta))
  {
    prec_alpha = alpha
    prec_beta  = beta
  }
  grid <- seq(0,max(s),length.out=lengthout+1)
  u <- diff(grid)
  field <- grid[-1]-u/2
  sgrid <- grid
  event_new <- 0
  time <- 0
  where <- 1
  E.factor <- 0
  for (j in 1:lengthout)
  {
    count <- sum(s>sgrid[j] & s<=sgrid[j+1]) 
    if (count>1)
    {
      points <- s[s>sgrid[j] & s<=sgrid[j+1]]
      u <- diff(c(sgrid[j],points))
      event_new <- c(event_new,event[(where):(where+count-1)])
      time <- c(time,rep(field[j],count))
      E.factor <- c(E.factor,coal.factor[where:(where+count-1)]*u)
      where <- where+count
      if (max(points)<sgrid[j+1])
      {
        event_new <- c(event_new,0)
        time <- c(time,field[j])
        E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-max(points)))
      }
    }
    if (count==1)
    {
      event_new <- c(event_new,event[where])
      points <- s[s>sgrid[j] & s<=sgrid[j+1]]
      if (points==sgrid[j+1])
      {
        E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
        time <- c(time,field[j])
        where <- where+1
      }
      else
      {
        event_new <- c(event_new,0)
        E.factor <- c(E.factor,coal.factor[where]*(points-sgrid[j]))
        E.factor <- c(E.factor,coal.factor[where+1]*(sgrid[j+1]-points))
        time <- c(time,rep(field[j],2))
        where <- where+1
      }
    }
    if (count==0)
    {
      event_new <- c(event_new,0)
      E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
      time <- c(time,field[j])
    }
    
  }
  time2 <- time
  event_new2 <- event_new
  E.factor2 <- E.factor
  
  for (j in 1:lengthout)
  {
    count <- sum(time2==field[j])
    if (count>1)
    {
      indic <- seq(1:length(event_new2))[time2==field[j]]
      if (sum(event_new2[indic])==0)
      {
        event_new2 <- event_new2[-indic[-1]]
        time2 <- time2[-indic[-1]]
        temp <- sum(E.factor2[indic])
        E.factor2[indic[1]] <- temp
        E.factor2 <- E.factor2[-indic[-1]]
      }
      #else {}
    }
  }
  
  #E.factor2[E.factor2 == 0] = exp(-1e6)
  E.factor2.log = log(E.factor2)
  E.factor2.log[E.factor2 == 0] = log_zero
  #print(E.factor2.log)
  
  #   E.factor2[1:34] <- rep(1,34)
  data <- list(y=event_new2[-1],event=event_new2[-1],time=time2[-1],E=E.factor2.log[-1])
  formula <- y~-1+f(time,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
  mod4 <- INLA::inla(formula,family="poisson",data=data,offset=E,control.predictor=list(compute=TRUE))
  
  return(list(result=mod4,grid=grid,data=data,E=E.factor2.log, x = field))
}

calculate_pref_old = function(coal.factor,s,event,lengthout,prec_alpha=0.01,prec_beta=0.01)
{
  grid <- seq(0,max(s),length.out=lengthout+1)
  u <- diff(grid)
  field <- grid[-1]-u/2
  sgrid <- grid
  event_new <- 0
  time <- 0
  where <- 1
  E.factor <- 0
  for (j in 1:lengthout)
  {
    count <- sum(s>sgrid[j] & s<=sgrid[j+1]) 
    if (count>1)
    {
      points <- s[s>sgrid[j] & s<=sgrid[j+1]]
      u <- diff(c(sgrid[j],points))
      event_new <- c(event_new,event[(where):(where+count-1)])
      time <- c(time,rep(field[j],count))
      E.factor <- c(E.factor,coal.factor[where:(where+count-1)]*u)
      where <- where+count
      if (max(points)<sgrid[j+1])
      {
        event_new <- c(event_new,0)
        time <- c(time,field[j])
        E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-max(points)))
      }
    }
    if (count==1)
    {
      event_new <- c(event_new,event[where])
      points <- s[s>sgrid[j] & s<=sgrid[j+1]]
      if (points==sgrid[j+1])
      {
        E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
        time <- c(time,field[j])
        where <- where+1
      }
      else
      {
        event_new <- c(event_new,0)
        E.factor <- c(E.factor,coal.factor[where]*(points-sgrid[j]))
        E.factor <- c(E.factor,coal.factor[where+1]*(sgrid[j+1]-points))
        time <- c(time,rep(field[j],2))
        where <- where+1
      }
    }
    if (count==0)
    {
      event_new <- c(event_new,0)
      E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
      time <- c(time,field[j])
    }
    
  }
  time2 <- time
  event_new2 <- event_new
  E.factor2 <- E.factor
  
  for (j in 1:lengthout)
  {
    count <- sum(time2==field[j])
    if (count>1)
    {
      indic <- seq(1:length(event_new2))[time2==field[j]]
      if (sum(event_new2[indic])==0)
      {
        event_new2 <- event_new2[-indic[-1]]
        time2 <- time2[-indic[-1]]
        temp <- sum(E.factor2[indic])
        E.factor2[indic[1]] <- temp
        E.factor2 <- E.factor2[-indic[-1]]
      }
      #else {}
    }
  }
  
  #   E.factor2[1:34] <- rep(1,34)
  #data <- list(y=event_new2[-1],event=event_new2[-1],time=time2[-1],E=log(E.factor2[-1]))
  #formula <- y~-1+f(time,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
  #mod4 <- inla(formula,family="poisson",data=data,offset=E,control.predictor=list(compute=TRUE),...)
  
  n1 <- length(event_new2[-1])
  n2 <- length(field)
  Y  <- matrix(NA,n1+n2,2)
  
  #dd <- heterochronous.gp.stat(influenza.tree)
  #s.time <- dd$sample.times
  #n.sample <- dd$sampled.lineages
  newcount <- rep(0,length(grid)-1)
  
  # MK: added to replace reading from influenza.tree above
  samps = s[-1][event==0]
  
  for (j in 1:length(newcount))
  {
    # MK: altered to samps version
    newcount[j] <- sum(samps > grid[j] & samps <= grid[j+1])
  }
  newcount[1] <- newcount[1]+1
  newcount[(sum(grid<max(samps))+1):lengthout] <- NA
  
  data.sampling2<-data.frame(y=newcount,time=field,E=log(diff(grid)))
  formula.sampling2=y~1+f(time,model="rw1",hyper = list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
  
  # MK: Added during functionizing
  E = diff(s) * coal.factor
  E[1]=1
  
  mod.sampling2<-inla(formula.sampling2,family="poisson",data=data.sampling2,offset=E,control.predictor=list(compute=TRUE))
  
  return(list(result=mod.sampling2,grid=grid))
}

calculate_moller_hetero_pref_old <- function(coal.factor,s,event,lengthout,prec_alpha=0.01,prec_beta=0.01,beta1_prec = 0.001,log_zero=-100,alpha=NULL,beta=NULL)
{
  if (prec_alpha == 0.01 & prec_beta==0.01 & !is.null(alpha) & !is.null(beta))
  {
    prec_alpha = alpha
    prec_beta  = beta
  }
  
  grid <- seq(0,max(s),length.out=lengthout+1)
  u <- diff(grid)
  field <- grid[-1]-u/2
  sgrid <- grid
  event_new <- 0
  time <- 0
  where <- 1
  E.factor <- 0
  for (j in 1:lengthout)
  {
    count <- sum(s>sgrid[j] & s<=sgrid[j+1]) 
    if (count>1)
    {
      points <- s[s>sgrid[j] & s<=sgrid[j+1]]
      u <- diff(c(sgrid[j],points))
      event_new <- c(event_new,event[(where):(where+count-1)])
      time <- c(time,rep(field[j],count))
      E.factor <- c(E.factor,coal.factor[where:(where+count-1)]*u)
      where <- where+count
      if (max(points)<sgrid[j+1])
      {
        event_new <- c(event_new,0)
        time <- c(time,field[j])
        E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-max(points)))
      }
    }
    if (count==1)
    {
      event_new <- c(event_new,event[where])
      points <- s[s>sgrid[j] & s<=sgrid[j+1]]
      if (points==sgrid[j+1])
      {
        E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
        time <- c(time,field[j])
        where <- where+1
      }
      else
      {
        event_new <- c(event_new,0)
        E.factor <- c(E.factor,coal.factor[where]*(points-sgrid[j]))
        E.factor <- c(E.factor,coal.factor[where+1]*(sgrid[j+1]-points))
        time <- c(time,rep(field[j],2))
        where <- where+1
      }
    }
    if (count==0)
    {
      event_new <- c(event_new,0)
      E.factor <- c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
      time <- c(time,field[j])
    }
    
  }
  time2 <- time
  event_new2 <- event_new
  E.factor2 <- E.factor
  
  for (j in 1:lengthout)
  {
    count <- sum(time2==field[j])
    if (count>1)
    {
      indic <- seq(1:length(event_new2))[time2==field[j]]
      if (sum(event_new2[indic])==0)
      {
        event_new2 <- event_new2[-indic[-1]]
        time2 <- time2[-indic[-1]]
        temp <- sum(E.factor2[indic])
        E.factor2[indic[1]] <- temp
        E.factor2 <- E.factor2[-indic[-1]]
      }
      #else {}
    }
  }
  
  #   E.factor2[1:34] <- rep(1,34)
  #data <- list(y=event_new2[-1],event=event_new2[-1],time=time2[-1],E=log(E.factor2[-1]))
  #formula <- y~-1+f(time,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
  #mod4 <- inla(formula,family="poisson",data=data,offset=E,control.predictor=list(compute=TRUE),...)
  
  n1 <- length(event_new2[-1])
  n2 <- length(field)
  Y  <- matrix(NA,n1+n2,2)
  
  #dd <- heterochronous.gp.stat(influenza.tree)
  #s.time <- dd$sample.times
  #n.sample <- dd$sampled.lineages
  newcount <- rep(0,length(grid)-1)
  
  # MK: added to replace reading from influenza.tree above
  samps = s[-1][event==0]
  
  for (j in 1:length(newcount))
  {
    # MK: altered to samps version
    newcount[j] <- sum(samps > grid[j] & samps <= grid[j+1])
  }
  newcount[1] <- newcount[1]+1
  newcount[(sum(grid<max(samps))+1):lengthout] <- NA
  
  #data.sampling2<-data.frame(y=newcount,time=field,E=log(diff(grid)))
  #formula.sampling2=y~1+f(time,model="rw1",hyper = list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)
  
  # MK: Added during functionizing
  E = diff(s) * coal.factor
  E[1]=1
  
  #print("Got here 1")
  
  beta0<-c(rep(0,n1),rep(1,n2))
  newE<-rep(NA,n1+n2)
  
  #MK: added to resolve INLA failure
  E.factor2.log = log(E.factor2)
  E.factor2.log[E.factor2 == 0] = log_zero
  newE[1:n1]<-E.factor2.log[-1]
  #newE[1:n1]<-log(E.factor2[-1])
  
  newE[(n1+1):(n1+n2)]<-log(diff(grid))
  megafield<-c(time2[-1],field)
  
  #print("Got here 2")
  
  Y[1:n1,1]<-event_new2[-1]
  #print("Got here 3")
  Y[(n1+1):(n2+n1),2]<-newcount
  #print("Got here 4")
  r<-c(rep(1,n1),rep(2,n2))
  w<-c(rep(1,n1),rep(-1,n2))
  # data.pref<-data.frame(Y=Y,beta0=beta0,r=r,megafield=megafield,E=newE,w=w)
  # formula.pref.rep<-Y~-1+beta0+f(megafield,w,model="rw1",replicate=r,hyper=list(prec = list(param = c(.001, .001))),constr=FALSE)
  ii<-c(time2[-1],rep(NA,n2))
  jj<-c(rep(NA,n1),field)
  
  #print("Got here 5")
  
  formula.pref<-Y~-1+beta0+f(ii,model="rw1",hyper=list(prec = list(param = c(prec_alpha, prec_beta))),constr=FALSE)+
    f(jj,w,copy="ii",fixed=FALSE,param=c(0,beta1_prec))
  #print("Got here 6")
  
  # MK: Changed data.frame to list
  data.pref<-list(Y=Y,beta0=beta0,ii=ii,jj=jj,E=newE,w=w)
  
  #print("Got here 7")
  #print(data.pref)
  mod.pref<-INLA::inla(formula.pref,family=c("poisson","poisson"),offset=E,data=data.pref,control.predictor=list(compute=TRUE))
  
  #print("Got here 8")
  
  return(list(result=mod.pref,grid=grid, x = field))
}
