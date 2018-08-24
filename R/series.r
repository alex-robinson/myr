# To do: merge functions in here that are elsewhere (derivative, etc)

#' @export
runmean = function(x,y,dt)
{ # Use function runmed to calculate running mean over a window of width dt, 
  # but ensure that time steps are evenly spaced.

  dt_min = min(diff(x))
  x1 = seq(min(x)-10*dt_min,max(x)+10*dt_min,by=dt_min)
  y1 = approx(x,y,x1,rule=2)$y 

  k = max(1,round(dt/dt_min)) %/% 2 *2 + 1

  y1sm = runmed(y1,k=k)

  ysm  = approx(x1,y1sm,x)$y

  return(ysm)
}

#' @export
my.loess = function(x,y,L)
{ # Calculate smoothing for a given window half-width L
  # (assumes even spacing in x)

  npts = (L*2+1)/diff(x)[1]
  ysm = predict(loess(y~x,span=npts/length(x)))

  return(ysm)
}

