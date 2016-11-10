float2gray = function(f)
{
  return(sprintf("#%02x%02x%02x", floor((1-f) * 255), floor((1-f) * 255), floor((1-f) * 255)))
}

#' Plot heatmap of a histogram
#' 
#' @param hist \code{histogram} object to be displayed.
#'
#' @param y numeric y-coordinate to display heatmap.
#' @param wd numeric width of heatmap in y-units.
#'
#' @export
hist2heat = function(hist, y, wd)
{
  breaks = hist$breaks
  counts = hist$counts
  upper  = max(counts)
  n = length(counts)
  cols = float2gray(counts / upper)
  graphics::segments(x0 = breaks[1:n], y0=y, x1 = breaks[2:(n+1)], y1=y, lwd=wd, col=cols, lend=1)
}
