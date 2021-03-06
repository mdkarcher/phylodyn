% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simple_models.R
\name{ml_exp_coal}
\alias{ml_exp_coal}
\title{Maximum likelihood exponential model}
\usage{
ml_exp_coal(samp_times, n_sampled, coal_times, betas = c(0, 0))
}
\arguments{
\item{samp_times}{vector of sampling times.}

\item{n_sampled}{vector of number sampled at \code{samp_times}.}

\item{coal_times}{vector of coalescent times.}

\item{betas}{starting values for optimization.}
}
\description{
Maximum likelihood exponential model
}
