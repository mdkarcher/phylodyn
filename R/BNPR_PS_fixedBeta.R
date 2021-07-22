infer_coal_samp_fixed <- function (samp_times, coal_times, n_sampled = NULL, fns = NULL, 
          lengthout = 100, prec_alpha = 0.01, prec_beta = 0.01, beta1_prec = 0.001, 
          use_samp = FALSE, log_fns = TRUE, simplify = FALSE, events_only = FALSE, 
          derivative = FALSE, link = 1) 
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("INLA needed for this function to work. Use install.packages(\"INLA\", repos=c(getOption(\"repos\"), INLA=\"https://inla.r-inla-download.org/R/stable\"), dep=TRUE).", 
         call. = FALSE)
  }
  if (min(coal_times) < min(samp_times)) 
    stop("First coalescent time occurs before first sampling time")
  if (max(samp_times) > max(coal_times)) 
    stop("Last sampling time occurs after last coalescent time")
  grid <- seq(min(samp_times), max(coal_times), length.out = lengthout + 
                1)
  if (is.null(n_sampled)) 
    n_sampled <- rep(1, length(samp_times))
  coal_data <- coal_stats(grid = grid, samp_times = samp_times, 
                          n_sampled = n_sampled, coal_times = coal_times)
  if (simplify) 
    coal_data <- with(coal_data, condense_stats(time = time, 
                                                event = event, E = E))
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  if (!use_samp) {
    data <- with(coal_data, data.frame(y = event, time = time, 
                                       E_log = E_log))
    formula <- y ~ -1 + f(time, model = "rw1", hyper = hyper, 
                          constr = FALSE)
    family <- "poisson"
  }
  else if (use_samp) {
    if (events_only) 
      samp_data <- samp_stats(grid = grid, samp_times = samp_times)
    else samp_data <- samp_stats(grid = grid, samp_times = samp_times, 
                                 n_sampled = n_sampled)
    data <- joint_stats(coal_data = coal_data, samp_data = samp_data)
    if (is.null(fns)) {
      formula <- Y ~ -1 + beta0 + f(time, model = "rw1", hyper = hyper, constr = FALSE) # + f(time2, w, copy = "time", fixed = FALSE, param = c(0, beta1_prec))
    }
    else {
      vals <- NULL
      bins <- sum(data$beta0 == 0)
      for (fni in fns) {
        if (log_fns) 
          vals <- cbind(vals, c(rep(0, bins), log(fni(samp_data$time))))
        else vals <- cbind(vals, c(rep(0, bins), fni(samp_data$time)))
      }
      data$fn <- vals
      formula <- Y ~ -1 + beta0 + fn + f(time, model = "rw1", hyper = hyper, constr = FALSE) # + f(time2, w, copy = "time", fixed = FALSE, param = c(0, beta1_prec))
    }
    family <- c("poisson", "poisson")
  }
  else stop("Invalid use_samp value, should be boolean.")
  if (derivative) {
    Imat <- diag(lengthout)
    A <- utils::head(Imat, -1) - utils::tail(Imat, -1)
    field <- grid[-1] - diff(grid)/2
    A <- diag(1/diff(field)) %*% A
    A[A == 0] <- NA
    lc_many <- INLA::inla.make.lincombs(time = A)
  }
  else {
    lc_many <- NULL
  }
  mod <- INLA::inla(formula, family = family, data = data, 
                    lincomb = lc_many, offset = data$E_log,
                    control.predictor = list(compute = TRUE, link = link),
                    control.inla = list(lincomb.derived.only = FALSE),
                    control.fixed = list(
                      mean = list(fn = 0),
                      prec = list(fn = 1e20) )
                     )
  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}
#environment(infer_coal_samp_fixed) <- environment(phylodyn:::infer_coal_samp)
#attributes(infer_coal_samp_fixed) <- attributes(phylodyn:::infer_coal_samp)

###
BNPR_fixed <- function (data, lengthout = 100, pref = TRUE, prec_alpha = 0.01, 
          prec_beta = 0.01, beta1_prec = 0.001, fns = NULL, log_fns = TRUE, 
          simplify = TRUE, derivative = FALSE, forward = TRUE, link = 1) 
{
  if (class(data) == "phylo") {
    phy <- summarize_phylo(data)
  }
  else if (all(c("coal_times", "samp_times", "n_sampled") %in% 
               names(data))) {
    phy <- with(data, list(samp_times = samp_times, coal_times = coal_times, 
                           n_sampled = n_sampled))
  }
  result <- infer_coal_samp_fixed(samp_times = phy$samp_times, coal_times = phy$coal_times, 
                            n_sampled = phy$n_sampled, fns = fns, lengthout = lengthout, 
                            prec_alpha = prec_alpha, prec_beta = prec_beta, beta1_prec = beta1_prec, 
                            use_samp = pref, log_fns = log_fns, simplify = simplify, 
                            derivative = derivative, link = link)
  result$samp_times <- phy$samp_times
  result$n_sampled <- phy$n_sampled
  result$coal_times <- phy$coal_times
  result$effpop <- exp(-result$result$summary.random$time$`0.5quant`)
  result$effpopmean <- exp(-result$result$summary.random$time$mean)
  result$effpop975 <- exp(-result$result$summary.random$time$`0.025quant`)
  result$effpop025 <- exp(-result$result$summary.random$time$`0.975quant`)
  result$summary <- with(result$result$summary.random$time, 
                         data.frame(time = ID, mean = exp(-mean), sd = sd * exp(-mean), 
                                    quant0.025 = exp(-`0.975quant`), quant0.5 = exp(-`0.5quant`), 
                                    quant0.975 = exp(-`0.025quant`)))
  if (derivative) {
    if (forward) 
      ind <- c(1:(lengthout - 1), (lengthout - 1))
    else ind <- c(1, 1:(lengthout - 1))
    result$derivative <- with(result$result$summary.lincomb, 
                              data.frame(time = result$x, mean = -mean[ind], sd = sd[ind], 
                                         quant0.025 = -`0.975quant`[ind], quant0.5 = -`0.5quant`[ind], 
                                         quant0.975 = -`0.025quant`[ind]))
  }
  if (pref) {
    result$beta0 <- result$result$summary.fixed["beta0", 
                                                "0.5quant"]
    result$beta0summ <- result$result$summary.fixed["beta0", 
                                                    ]
    rownames(result$beta0summ) <- "Beta0"
    result$beta1 <- result$result$summary.hyperpar[2, "0.5quant"]
    result$beta1summ <- result$result$summary.hyperpar[2, ]
    rownames(result$beta1summ) <- "Beta1"
  }
  return(result)
}
#environment(BNPR_fixed) <- environment(phylodyn::BNPR)
#attributes(BNPR_fixed) <- attributes(phylodyn::BNPR)