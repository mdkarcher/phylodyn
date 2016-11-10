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
  return(INLA::inla.models(...))
}
