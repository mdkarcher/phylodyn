shade_band = function(x, ylo, yhi, col="gray")
{
  polygon(c(x, rev(x)), c(yhi, rev(ylo)), col=col, border=NA)
}
