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
    stop('INLA needed for this function to work. Use install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable").',
         call. = FALSE)
  }
  
  return(INLA::inla.models(...))
}
