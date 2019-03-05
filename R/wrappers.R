#' INLA wrapper
#' 
#' Wrapper for \code{INLA::inla.models} to work around a potential namespace bug
#' in \code{INLA::inla}.
#' 
#' @param ... passes all arguments to \code{INLA::inla.models}.
#'   
#' @export
inla.models <- function(...)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE).',
         call. = FALSE)
  }
  
  return(INLA::inla.models(...))
}
