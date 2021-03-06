% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_check.R
\name{samp_check}
\alias{samp_check}
\title{Predictive check for sampling model}
\usage{
samp_check(samp_times, n_sampled, logpop, grid, betas, covariates = NULL,
  cov_betas = NULL, pow_covariates = NULL, pow_cov_betas = NULL)
}
\arguments{
\item{samp_times}{numeric vector of sampling times.}

\item{n_sampled}{integer vector of lineages sampled at each sampling time.}

\item{logpop}{numeric vector of log-effective population sizes.}

\item{grid}{numeric grid of time points.}

\item{betas}{numeric vector of log-linear coefficients for sampling model.}

\item{covariates}{list of vectors representing covariate values.}

\item{cov_betas}{numeric coefficients for each element of covariates.}

\item{pow_covariates}{list of vectors representing interaction covariate
values.}

\item{pow_cov_betas}{numeric coefficients for each element of pow_covariates.}
}
\value{
named numeric vector: "obs" contains the chi-square discrepancy for
  the supplied sampling times, "rep" contains å chi-square discrepancy for
  the replicated sampling times.
}
\description{
Predictive check for sampling model
}
