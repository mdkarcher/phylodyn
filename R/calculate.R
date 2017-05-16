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
  lineages <- utils::head(cumsum(lineage_change), -1) # remove entry for the post-final-coalescent-event open interval
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
#'   times \code{coal_times}, sampling times \code{samp_times}, and number 
#'   sampled per sampling time \code{n_sampled}.
#' @param lengthout numeric specifying number of grid points.
#' @param pref logical. Should the preferential sampling model be used?
#' @param prec_alpha,prec_beta numerics specifying gamma prior for precision 
#'   \eqn{\tau}.
#' @param beta1_prec numeric specifying precision for normal prior on 
#'   \eqn{\beta_1}.
#' @param fns list containing functions of covariates.
#' @param log_fns logical whether or not to to apply a log-transformation to
#'   the output of the functions in \code{fns}.
#' @param simplify logical whether to fully bucket all Poisson points.
#' @param derivative logical whether to calculate estimates of the 
#'   log-derivative.
#' @param forward logical whether to use the finite difference approximations of
#'   the log-derivative as a forward or backward derivative.
#'   
#' @return Phylodynamic reconstruction of effective population size at grid 
#'   points. \code{result} contains the INLA output, \code{data} contains the 
#'   information passed to INLA, \code{grid} contains the grid end points, 
#'   \code{x} contains the grid point centers, \code{effpop} contains a vector 
#'   of the posterior median effective population size estimates, 
#'   \code{effpop025} and \code{effpop975} contain the 2.5th and 97.5th 
#'   posterior percentiles, \code{summary} contains a data.frame of the 
#'   estimates, and \code{derivative} (if \code{derivative = TRUE}) contains a
#'   data.frame summarizing the log-derivative.
#' @export
#' 
#' @examples
#' data("NY_flu")
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'  res = BNPR(NY_flu)
#'  plot_BNPR(res)
#' }
BNPR <- function(data, lengthout = 100, pref=FALSE, prec_alpha=0.01,
                 prec_beta=0.01, beta1_prec = 0.001, fns = NULL, log_fns = TRUE,
                 simplify = TRUE, derivative = FALSE, forward = TRUE)
{
  if (class(data) == "phylo")
  {
    phy <- summarize_phylo(data)
  }
  else if (all(c("coal_times", "samp_times", "n_sampled") %in% names(data)))
  {
    phy <- with(data, list(samp_times = samp_times, coal_times = coal_times,
                           n_sampled = n_sampled))
  }
  
  result <- infer_coal_samp(samp_times = phy$samp_times, coal_times = phy$coal_times,
                            n_sampled = phy$n_sampled, fns = fns, lengthout = lengthout,
                            prec_alpha = prec_alpha, prec_beta = prec_beta,
                            beta1_prec = beta1_prec, use_samp = pref, log_fns = log_fns,
                            simplify = simplify, derivative = derivative)
  
  result$samp_times <- phy$samp_times
  result$n_sampled  <- phy$n_sampled
  result$coal_times <- phy$coal_times
  
  result$effpop     <- exp(-result$result$summary.random$time$`0.5quant`)
  result$effpopmean <- exp(-result$result$summary.random$time$mean)
  result$effpop975  <- exp(-result$result$summary.random$time$`0.025quant`)
  result$effpop025  <- exp(-result$result$summary.random$time$`0.975quant`)
  
  result$summary <- with(result$result$summary.random$time,
                         data.frame(time = ID, mean = exp(-mean),
                                    sd = sd * exp(-mean),
                                    quant0.025 = exp(-`0.975quant`),
                                    quant0.5 = exp(-`0.5quant`),
                                    quant0.975 = exp(-`0.025quant`)))
  
  if (derivative)
  {
    if (forward)
      ind <- c(1:(lengthout-1), (lengthout-1))
    else
      ind <- c(1, 1:(lengthout-1))
    
    result$derivative <- with(result$result$summary.lincomb,
                              data.frame(time = result$x, mean = -mean[ind], sd = sd[ind],
                                         quant0.025 = -`0.975quant`[ind],
                                         quant0.5   = -`0.5quant`[ind],
                                         quant0.975 = -`0.025quant`[ind]))
  }
  
  if (pref)
  {
    result$beta0     <- result$result$summary.fixed["beta0","0.5quant"]
    result$beta0summ <- result$result$summary.fixed["beta0",]
    rownames(result$beta0summ) <- "Beta0"
    result$beta1     <- result$result$summary.hyperpar[2,"0.5quant"]
    result$beta1summ <- result$result$summary.hyperpar[2,]
    rownames(result$beta1summ) <- "Beta1"
  }
  
  return(result)
}

#' @describeIn BNPR Uses preferential sampling model.
#' @export
BNPR_PS <- function(data, lengthout = 100, prec_alpha=0.01, prec_beta=0.01,
                    beta1_prec = 0.001, fns = NULL, log_fns = TRUE,
                    simplify = TRUE, derivative = FALSE, forward = TRUE)
{
  return(BNPR(data = data, lengthout = lengthout, pref = TRUE,
              prec_alpha = prec_alpha, prec_beta = prec_beta,
              beta1_prec = beta1_prec, fns = fns, log_fns = log_fns,
              simplify = simplify, derivative = derivative, forward = forward))
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
  
  Cfun <- stats::stepfun(x = s, y = c(0, coal_factor, 0), right = TRUE)
  Cvec <- Cfun(sgrid[-1])
  E <- diff(sgrid)*Cvec
  
  E_log = log(E)
  E_log[E == 0] = log_zero
  
  return(data.frame(time = time, event = event_out[-1], E = E, E_log = E_log))
}

condense_stats <- function(time, event, E, log_zero = -100)
{
  result <- stats::aggregate(event ~ time, FUN = sum)
  result$E <- stats::aggregate(E ~ time, FUN = sum)$E
  
  E_log = log(result$E)
  E_log[result$E == 0] = log_zero
  result$E_log <- E_log
  
  return(result)
}

infer_coal <- function(samp_times, coal_times, n_sampled = NULL, lengthout = 100,
                       prec_alpha = 0.01, prec_beta = 0.01, simplify = FALSE,
                       derivative = FALSE)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable").',
         call. = FALSE)
  }
  
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
    coal_data <- with(coal_data, condense_stats(time = time, event = event, E=E))
  
  data <- with(coal_data, data.frame(y = event, time = time, E_log = E_log))
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  formula <- y ~ -1 + f(time, model="rw1", hyper = hyper, constr = FALSE)
  
  if (derivative)
  {
    Imat <- diag(lengthout)
    A <- utils::head(Imat, -1) - utils::tail(Imat, -1)
    field <- grid[-1] - diff(grid)/2
    A <- diag(1/diff(field)) %*% A
    A[A == 0] <- NA
    
    lc_many <- INLA::inla.make.lincombs(time = A)
    
    mod <- INLA::inla(formula, family = "poisson", data = data, lincomb = lc_many,
                      control.predictor = list(compute=TRUE),
                      control.inla = list(lincomb.derived.only=FALSE))
  }
  else
  {
    mod <- INLA::inla(formula, family = "poisson", data = data, offset = data$E_log,
                      control.predictor = list(compute=TRUE))
  }
  
  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}

samp_stats <- function(grid, samp_times, n_sampled = NULL, trim_end = FALSE)
{
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2
  E <- diff(grid)
  
  bins <- cut(x = samp_times, breaks = grid, include.lowest = TRUE)
  
  if (is.null(n_sampled))
    count <- as.vector(table(bins))
  else
  {
    tab <- stats::aggregate(n_sampled ~ bins, FUN = sum, labels = FALSE)
    count <- rep(0, lengthout)
    count[as.numeric(tab$bins)] <- tab$n_sampled
  }
  
  count[utils::head(grid, -1) >= max(samp_times)] <- NA
  result <- data.frame(time = field, count = count, E = E, E_log = log(E))
  
  if (trim_end)
    result <- result[stats::complete.cases(result),]
  
  return(result)
}

infer_samp <- function(samp_times, n_sampled = NULL, lengthout = 100,
                       prec_alpha = 0.01, prec_beta = 0.01)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable").',
         call. = FALSE)
  }
  
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

joint_stats <- function(coal_data, samp_data)
{
  n1 <- length(coal_data$time)
  n2 <- length(samp_data$time)
  beta0 <- c(rep(0, n1), rep(1, n2))
  E_log <- c(coal_data$E_log, samp_data$E_log)
  Y <- matrix(c(coal_data$event, rep(NA, n2), rep(NA, n1), samp_data$count),
              nrow = n1 + n2, byrow = FALSE)
  w <- c(rep(1, n1), rep(-1, n2))
  time  <- c(coal_data$time, rep(NA, n2))
  time2 <- c(rep(NA, n1), samp_data$time)
  
  return(list(Y = Y, beta0 = beta0, time = time, time2 = time2, w = w, E_log = E_log))
}

infer_coal_samp <- function(samp_times, coal_times, n_sampled=NULL, fns = NULL,
                            lengthout=100, prec_alpha=0.01, prec_beta=0.01,
                            beta1_prec=0.001, use_samp = FALSE, log_fns = TRUE,
                            simplify = FALSE, events_only = FALSE,
                            derivative = FALSE)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable").',
         call. = FALSE)
  }
  
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
  
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  
  if (!use_samp)
  {
    data <- with(coal_data, data.frame(y = event, time = time, E_log = E_log))
    
    formula <- y ~ -1 + f(time, model="rw1", hyper = hyper, constr = FALSE)
    family <- "poisson"
  }
  else if (use_samp)
  {
    if (events_only)
      samp_data <- samp_stats(grid = grid, samp_times = samp_times)
    else
      samp_data <- samp_stats(grid = grid, samp_times = samp_times,
                              n_sampled = n_sampled)
    
    data <- joint_stats(coal_data = coal_data, samp_data = samp_data)
    
    if (is.null(fns))
    {
      formula <- Y ~ -1 + beta0 +
        f(time, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))
    }
    else
    {
      vals <- NULL
      bins <- sum(data$beta0 == 0)
      for (fni in fns)
      {
        if (log_fns)
          vals <- cbind(vals, c(rep(0, bins), log(fni(samp_data$time))))
        else
          vals <- cbind(vals, c(rep(0, bins), fni(samp_data$time)))
      }
      data$fn <- vals
      
      formula <- Y ~ -1 + beta0 + fn +
        f(time, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))
    }
    
    family <- c("poisson", "poisson")
  }
  else
    stop("Invalid use_samp value, should be boolean.")
  
  if (derivative)
  {
    Imat <- diag(lengthout)
    A <- utils::head(Imat, -1) - utils::tail(Imat, -1)
    field <- grid[-1] - diff(grid)/2
    A <- diag(1/diff(field)) %*% A
    A[A == 0] <- NA
    
    lc_many <- INLA::inla.make.lincombs(time = A)
  }
  else
  {
    lc_many <- NULL
  }
  
  mod <- INLA::inla(formula, family = family, data = data,
                    lincomb = lc_many, offset = data$E_log,
                    control.predictor = list(compute=TRUE),
                    control.inla = list(lincomb.derived.only=FALSE))
  
  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}

infer_coal_deriv <- function(samp_times, coal_times, n_sampled = NULL, lengthout = 100,
                             prec_alpha = 0.01, prec_beta = 0.01, simplify = FALSE)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable").',
         call. = FALSE)
  }
  
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
  formula <- y ~ -1 + f(time, model="rw1", hyper = hyper, constr = FALSE) + offset(data$E_log)
  
  Imat <- diag(lengthout)
  A <- utils::head(Imat, -1) - utils::tail(Imat, -1)
  A[A == 0] <- NA
  
  #lcmat <- rbind(c(1, -1, rep(NA, lengthout - 2)), c(NA, 1, -1, rep(NA, lengthout - 3)))
  lc_many <- INLA::inla.make.lincombs(time = A)
  #names(lc1) = "lc1"
  
  mod <- INLA::inla(formula, family = "poisson", data = data, lincomb = lc_many,
                    control.predictor = list(compute=TRUE),
                    control.inla = list(lincomb.derived.only=FALSE))
  
  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}

infer_samp_exper <- function(samp_times, fns, n_sampled = NULL, lengthout = 100,
                             prec_alpha = 0.01, prec_beta = 0.01)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable").',
         call. = FALSE)
  }
  
  grid <- seq(min(samp_times),max(samp_times),length.out=lengthout+1)
  
  samp_data <- samp_stats(grid = grid, samp_times = samp_times,
                          n_sampled = n_sampled)
  
  vals = NULL
  for (fni in fns)
  {
    vals = cbind(vals, log(fni(samp_data$time)))
  }
  
  data <- with(samp_data, list(y = count, time = time, E_log = E_log))
  data$fn=vals
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  formula_sampling <- y ~ 1 + fn + f(time, model="rw1", hyper = hyper, constr=FALSE)
  
  mod <- INLA::inla(formula_sampling, family="poisson", data=data,
                    offset=data$E_log, control.predictor=list(compute=TRUE))
  
  return(list(result = mod, data = data, grid = grid, x = samp_data$time))
}

infer_coal_samp_exper <- function(samp_times, coal_times, n_sampled=NULL, fns = NULL,
                                  lengthout=100, prec_alpha=0.01, prec_beta=0.01,
                                  beta1_prec=0.001, use_samp = FALSE, log_fns = TRUE,
                                  simplify = FALSE, events_only = FALSE)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable").',
         call. = FALSE)
  }
  
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
  
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  
  if (events_only)
    samp_data <- samp_stats(grid = grid, samp_times = samp_times)
  else
    samp_data <- samp_stats(grid = grid, samp_times = samp_times,
                            n_sampled = n_sampled)
  
  data <- joint_stats(coal_data = coal_data, samp_data = samp_data)
  
  vals <- NULL
  bins <- sum(data$beta0 == 0)
  for (fni in fns)
  {
    vals <- cbind(vals, c(rep(0, bins), log(fni(samp_data$time))))
  }
  data$fn <- vals
  
  formula <- Y ~ -1 + beta0 + fn +
    f(time, model="rw1", hyper = hyper, constr = FALSE) +
    f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))
  
  family <- c("poisson", "poisson")
  
  mod <- INLA::inla(formula, family = family, data = data, offset = data$E_log,
                    control.predictor = list(compute=TRUE),
                    control.inla = list(lincomb.derived.only=FALSE))
  
  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}