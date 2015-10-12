# gen_Nehat = function(INLA_out, quant="0.5quant")
# {
#   out = INLA_out$result$summary.random$time
#   return(approxfun(x=out$ID, y=exp(-out[[quant]])))
# }
# 
# approx_INLA = function(INLA_out, xout, name="time")
# {
#   mod = INLA_out$result$summary.random[[name]]
#   x = mod$ID
#   y = exp(-mod$"0.5quant")
#   return(approx(x, y, xout, rule=2:1)$y)
# }
# 
# approx_INLA_lo = function(INLA_out, xout, name="time")
# {
#   mod = INLA_out$result$summary.random[[name]]
#   x = mod$ID
#   y = exp(-mod$"0.975quant")
#   return(approx(x, y, xout, rule=2:1)$y)
# }
# 
# approx_INLA_hi = function(INLA_out, xout, name="time")
# {
#   mod = INLA_out$result$summary.random[[name]]
#   x = mod$ID
#   y = exp(-mod$"0.025quant")
#   return(approx(x, y, xout, rule=2:1)$y)
# }
