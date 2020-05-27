
#' @export
rmse <- function(x,ii=c(1:length(x)),round=NA)
{
  # Filter for desired points
  x <- x[ii]
  
  ii1 <- which(!is.na(x))
  rmse <- sqrt(sum(x[ii1]^2)/length(x[ii1]))
  
  if (!is.na(round)) rmse <- round(rmse,round)
  
  return(rmse)
  
}

#' @export
rsquared <- function(x,xfit,ii=c(1:length(x)),round=NA)
{ # Calculate the R-squared value: R-squared = 1 - MSE/VAR(Y)

  # Filter for desired points
  ii0 <- c(1:length(x))
  ii1 <- which(!is.na(x) & !is.na(xfit) & ii0 %in% ii)
  x <- x[ii1]
  xfit <- xfit[ii1]
  
  xm <- mean(x)
  
  varx <- sum((x-xm)^2)
  
  MSE <- sum((x-xfit)^2)
  
  rsq <- 1 - MSE/varx
  
  if (!is.na(round)) rsq <- round(rsq,round)
  
  return(rsq)
  
}

#' @export
normalize <- function(x,sd1=NA,ave=NA)
{
  if (is.na(sd1)) sd1 <- sd(x,na.rm=TRUE)
  if (is.na(ave)) ave <- mean(x,na.rm=TRUE)
  
  x1 <- (x-ave)/sd1
  return(x1)
}

#' @export
standardize <- function(pres,ibase=c(1:length(pres)))
{ # Standardize the data of a given month
  
  # Determine the average pressure during the base time period
  # Then calculate the anomalies from this average pressure
  ave <- mean(pres[ibase],na.rm=TRUE)
  dpres <- pres - ave
  
  # Calculate the standard deviation of the anomalies
  sd1 <- sd(dpres[ibase],na.rm=TRUE)
  
  # Divide by the standard deviation of the anomalies to normalize (standardize) it
  stdpres <- dpres / sd1
  
  return(stdpres)
}

#' @export
derivative <- function(x,y,L=NULL)
{
    dy = c(NA,diff(y)/diff(x))
    dy[1] = dy[2]

    if (!is.null(L)) {  # Smooth the result using loess where L is window half-length
        npts = (L*2+1)/diff(x)[1]
        dy = predict(loess(dy~x,span=npts/length(x)))
    }

    return(dy)
}

check_max = function(xnow,x,max)
{   # Check if the nearest x is within the max allowed distance
    check = 1 
    if ( min(abs(x-xnow)>max,na.rm=TRUE) ) check = NA
    return(check)
}

#' @export
resample = function(x,y,xout,dtmax=NULL)
{   # Resample dataset, introduce NA values if necessary
    new = approx(x,y,xout=xout,rule=2)$y 
    if (!is.null(dtmax)) for (k in 1:length(xout)) new[k] = new[k]*check_max(xout[k],x,max=dtmax)
    return(new)
}

## Adding histogram to a plot
#' @export
my.hist <- function(hi,freq=TRUE,b=FALSE,col='grey95',border='grey70',lwd=1,lty=1,filled="standard")
{
  cat("myr:: Note that my.hist is deprecated. Please use myhist().","\n")

  x  <- hi$mids
  x2 <- hi$breaks
  y  <- hi$counts
  if (!freq) y <- hi$density
  
  if ( b ) { # USE NEW METHOD DATA
    x  <- hi$mids2
    x2 <- hi$breaks2
    y  <- hi$counts2
    if (!freq) y <- hi$density2
  }

  ## Standard histogram plot
  # for ( i in 1:length(x) ) {
  #   dx <- (x2[i+1]-x2[i]) / 2
  #   rect(x[i]-dx,0,x[i]+dx,y[i],col=col,border=border,lwd=lwd,lty=lty)
  # }
  
  ## histogram outline without vertical lines between bins
   n <- length(x2)
   xfill <- c(x2[1],x2[1])
   yfill <- c(0,y[1])
   for ( i in 1:n ) {
     xfill <- c(xfill,x2[i],x2[i+1],x2[i+1],x2[i+1])
     yfill <- c(yfill,y[i],y[i],y[i],y[i+1])
   }
   xfill <- c(xfill,x2[n],x2[n])
   yfill <- c(yfill,y[n],0)
   polygon(xfill,yfill,col=col,border=NA)
   lines(xfill,yfill,col=border,lwd=lwd,lty=lty)
  
  ## Polygon around values
  # polygon(hi$density.call,col=col,border=border,lwd=lwd,lty=lty)
  #polygon(c(x[1],x,x[length(x)]),c(0,y,0),col=col,border=border,lwd=lwd,lty=lty)

}

## Adding histogram to a plot bug fixes
#' @export
myhist = function (hi,freq=TRUE,col="grey95",border="grey70",lwd=1,lty=1,vertical.lines=TRUE) 
{
    x = hi$mids
    x2 = hi$breaks
    y = hi$counts
    if (!freq) y = hi$density
    n = length(x2)
    xfill = c(x2[1], x2[1])
    yfill = c(0, y[1])
    for (i in 2:(n-1)) {
        xfill = c(xfill,  x2[i], x2[i])
        yfill = c(yfill, y[i-1],  y[i])
    }
    xfill = c(xfill,  x2[n], x2[n])
    yfill = c(yfill, y[n-1], 0)

    polygon(xfill,yfill,col=col,border=NA)
    lines(xfill,yfill,col=border,lwd=lwd,lty=lty)

    if (vertical.lines) {
      # Include vertical lines separating bins
      segments(x0=x2[2:n],y0=0,x1=x2[2:n],y1=y,col=border,lwd=lwd*0.5,lty=lty)
    }
    
}

#' @export
get.conf <- function(x,weight,interval=95)
{ # Function to determine confidence/credence intervals
  # given a sample x and weights (weight should sum to 1)
  # Returns the indices that fall within the desired interval
  
  # Number of points
  n <- length(weight)

  # Order the weights from smallest to largest
  i1 <- order(weight)
  y1 <- weight[i1]

  # Loop over ordered weights, until interval is reached 
  # (ie, count until the interval of interest begins)
  w <- 0.0
  for ( i in 1:n) {
    w <- w + y1[i]
    if ( w >= (100-interval)*1e-2 ) break
  }
  
  # Collect the indices of points inside the desired interval
  # and the actual values
  ii95 <- i1[i:n]
  x95 <- x[ii95]

  # Return the indices
  return(ii95)
}

#' @export
get.threshold <- function(x,dens,percent)
{ # Calculate the level x1 at which sum(density) >= percent 

  x1  = rep(NA,length(percent))
  wgt = x1 

  # Number of points
  n <- length(dens)-1

  for (q in 1:length(x1)) {
    w = 0 
    for (i in 1:n) {
      w = w + (dens[i]*(x[i+1]-x[i]))
      if ( w >= percent[q]*1e-2 ) break 
    }
    x1[q] = x[i]
    wgt[q] = dens[i] 
  }

  return(data.frame(percent=percent,x=x1,weight=wgt))
}

#' @export
dnorm.integrated <- function(x,n=100,mean=0,sd=1)
{ # Get the integrated prob. density of bins (instead of just 1 value)
  # (function requested by Dr Rougier for statistics of hysteresis paper)

  x0 <- sort(unique(x))
  dd1 <- diff(x0)/2
  dd1 <- c(mean(dd1),dd1,mean(dd1))
  pp1 <- numeric(length(x0))
  p <- x*NA
  for ( q in 1:length(x0) ) {
    r0 <- x0[q]-dd1[q]
    r1 <- x0[q]+dd1[q+1]
    pp1[q] <- sum(dnorm(seq(r0,r1,length.out=n),mean=mean,sd=sd))/n
    p[which(x==x0[q])] <- pp1[q]
  }
  
  return(p)
}

#' @export
samples.syn <- function(x,p,n=1e3,dx=0.01,spar=0.6,from=range(x,na.rm=T)[1],to=range(x,na.rm=T)[2]) 
{ ## Method: Get spline of cdf and sample a large number of
  ## points from it. Return individual points with probabilities attached and
  ## and evenly spaced pdf as derivative of smoothed spline
  
  # Determine the sampled cumulative density
  ii   <- order(x)
  x1   <- x[ii]
  dx1  <- diff(x1)
  p1   <- p[ii]
  cdf1 <- cumsum(p1)
  
  # Add additional points to beginning and end of x to make spline 0:1
  xbuffer <- cumsum(rep(mean(dx1,na.rm=T)*10,10))
  x1 <- c(x1[1]-rev(xbuffer),x1,x1[length(x1)]+xbuffer)
  cdf1 <- c( rep(0.0,length(xbuffer)), cdf1, rep(1.0,length(xbuffer)) )
  
  # Calculate the spline function of the cdf
  cdf1fun <- splinefun(x=x1,y=cdf1,method="monoH.FC")
  
  ### FIRST calculate spline at random x points to get individual probabilities
  min1 = range(x[p>0])[1]; max1 = range(x[p>0])[2]
  xr <- sort(runif(n,min=min1,max=max1))
  cdfr <- cdf1fun(xr,deriv=0)
  cdfr[cdfr>1] <- 1
  cdfr[cdfr<0] <- 0

  # Calculate the probability of each point p[i] = cdf[i]-cdf[i-1]
  pr <- cdfr[2:length(cdfr)] - cdfr[1:(length(cdfr)-1)]
  pr <- c(0,pr)
  pr <- pr/sum(pr)
  
  ## NEXT calculate spline at evenly-spaced intervals to get the pdf and cdf
  mids <- seq(from=from,to=to,by=dx)
  cdf <- cdf1fun(mids,deriv=0)
  cdf[cdf>1] <- 1
  cdf[cdf<0] <- 0
  
  ## SMOOTH THE CDF
  ## Note: running mean doesn't work, because
  ## the cdf flattens out too much

  ## Smooth the spline now using smooth.spline
  ## (This wasn't possible originally, bc smooth.spline 
  ##  doesn't have option for monotonic splines)
  cdf = smooth.spline(x=mids,y=cdf,spar=spar)$y
  
  # Calculate the probability density
  # (Derivative of cdf, using central-differencing)
  pdf <- cdf*0
  for ( i in 3:(length(pdf)-2) ) {
    
    # Central difference (4th order error)
    pdf[i] <- ( -cdf[i+2] + 8*cdf[i+1] - 8*cdf[i-1] + cdf[i-2] ) / (12*diff(mids)[1])
  }
  pdf[pdf<0] = 0
  
  return(list(mids=mids,density=pdf,cdf=cdf,x1=x1,cdf1=cdf1,xr=xr,pr=pr))    

}


#' @export
my.density <- function(x,weights=rep(1,length(x)),subset=NULL,
                       from=NULL,to=NULL,n=512,dx=NA,dx.hist=0.1,
                       na.rm=T,type="kde",calchist=FALSE,verbose=TRUE)
{ # Function to calculate a PDF based on a sample of data x,
  # and given weights. Output should be similar to standard hist() function,
  # along with confidence intervals and a stats summary
  
  if (!is.null(subset)) {
    x       = x[subset]
    weights = weights[subset]
  }

  if (na.rm) {
    ii = which(!is.na(x+weights) )
    x       = x[ii]
    weights = weights[ii]
  }
  
  # Check range of x values
  lims0 = range(x)
  goodrange = TRUE
  if (diff(lims0)==0 | length(x) < 5) goodrange = FALSE

  if (is.null(from)) from = lims0[1]
  if (is.null(to))   to   = lims0[2]
  
  # If dx is given (desired output spacing), determine n
  if (!is.na(dx)) {
    tmp <- seq(from=from,to=to,by=dx)
    n   <- length(tmp)
  } else {
    tmp <- seq(from=from,to=to,length.out=n)
  }
  
  # Make sure weights are normalized
  weights = weights / sum(weights)
  
  if (!goodrange) { # Filler to allow function to keep going
    density.call = NA
    mids = tmp
    density = mids*0+1
    density = density / sum(density) / dx 
  } else if ( type=="kde" ) { # Use kernel density estimate
    density.call <- density(x,weights=weights,from=from,to=to,n=n)
    mids    <- density.call$x
    density <- density.call$y
  } else { # Use smoothed empirical CDF
    density.call <- samples.syn(x=x,p=weights,from=from,to=to,dx=dx)
    mids    <- density.call$mids
    density <- density.call$density
  }

  dx = diff(mids)[1]
  db = dx/2
  breaks  = seq(from=mids[1]-db,to=max(mids)+db,length.out=length(mids)+1)
  counts  = density*dx*length(x)
  
  ## NEW - calculate histogram instead of density
  hist = NA
  if (calchist) {
    hist = list()
    hist$breaks = seq(from=from,to=to,by=dx.hist)
    hist$mids   = hist$breaks[1:(length(breaks)-1)]+dx.hist/2
    hist$counts = hist$density = density*0
    for ( i in 1:length(hist$mids) ) {
      ii = which( x >= hist$breaks[i] & x < hist$breaks[i+1] )
      hist$density[i] = sum(weights[ii])
      hist$counts[i]  = length(ii)
    }
    hist$density = hist$density / diff(hist$breaks[1:2])
  
  }
  ## END NEW
  
  # Get various confidence intervals of interest (from density output)
  wts = density/sum(density); wts[is.na(wts)] = 0
  i95 = get.conf(mids,wts,interval=95)
  i90 = get.conf(mids,wts,interval=90)
  i66 = get.conf(mids,wts,interval=66)
  #i10 = get.conf(mids,wts,interval=10)
  #i01 = get.conf(mids,wts,interval=1)
  i10 = i01 = NA

  idens = get.conf(mids,wts,interval=5)
  best  = mean(mids[idens],na.rm=T)
  # best  = mids[which.max(density)]

  # Generate an output table (from density output)
  table =  data.frame(best=best,
                      r66a=range(mids[i66])[1],r66b=range(mids[i66])[2],
                      r90a=range(mids[i90])[1],r90b=range(mids[i90])[2],
                      r95a=range(mids[i95])[1],r95b=range(mids[i95])[2])
  
  table.new = c(table$r90a,table$r66a,best,table$r66b,table$r90b)

  if (!goodrange) table.new = rep(NA,5)
  
  # Get various confidence intervals of interest (from original data)
  # i95 = get.conf(x,weights,interval=95)
  # i90 = get.conf(x,weights,interval=90)
  # i66 = get.conf(x,weights,interval=66)
  # i10 = get.conf(x,weights,interval=10)
  # i01 = get.conf(x,weights,interval=1)
  # # Determine best (from top few points)
  # idens = get.conf(mids,density,interval=2)
  # best  = mean(mids[idens],na.rm=T)
  # # best  = mids[which.max(density)]

  # # Generate an output table (from original data)
  # table.orig = data.frame(best=best,
  #                         r66a=range(x[i66])[1],r66b=range(x[i66])[2],
  #                         #r90a=range(x[i90])[1],r90b=range(x[i90])[2],
  #                         r95a=range(x[i95])[1],r95b=range(x[i95])[2])
  #                         #r100a=range(x)[1],r100b=range(x)[2],
  #                         #bestsim=x[which.max(weights)])
  table.orig = NA 

  # Store all results
  out <- list(mids=mids,counts=counts,density=density,breaks=breaks,density.call=density.call,
              hist=hist,x=x,weights=weights,
              i95=i95,i90=i90,i66=i66,i10=i10,i01=i01,summary=table,summary.orig=table.orig,table=table.new)
  
  if (verbose) cat("Calculated my.density. Best x =",table$best,'\n')

  return(out)
}

#' @export
normalize_density <- function(dens,mids) 
{   # Returns the density renormalized (eg after multiplying probabilities together)

    dx = diff(mids[1:2])
    dens = dens / sum(dens) / dx
    return(dens)
}

#' @export
get_intervals <- function(dens,mids)
{   # Returns the 2.5, 33, 50, 66, 97.5 conf intervals

    wts = dens/sum(dens)
    wts[is.na(wts)] = 0
    i95 = get.conf(mids,wts,interval=95)
    i66 = get.conf(mids,wts,interval=66)

    idens = get.conf(mids,wts,interval=2)
    best  = mean(mids[idens],na.rm=T)

    # Generate an output table (from density output)
    table =  data.frame(best=best,
                        r66a=range(mids[i66])[1],r66b=range(mids[i66])[2],
                        r95a=range(mids[i95])[1],r95b=range(mids[i95])[2])
  
  table.new = c(table$r95a,table$r66a,best,table$r66b,table$r95b)

    return(table.new)
}

#' @export
join_pdfs <- function(dens2D,mids)
{   # Given density[mids,other], make a marginal pdf: density[mids]

    # Multiply all pdfs together
    n = dim(dens2D)[2]
    dens = dens2D[,1]
    for (k in 2:n) {
        dens = dens*dens2D[,k]
    }

    # Renormalize
    dens = normalize_density(dens,mids)

    out = list(x=mids,y=dens)
    return(out)
}

#' @export
join_pdfs1 <- function(dens2D,mids)
{   # Given density[mids,other], make a marginal pdf: density[mids]

    # Sum all pdfs together
    dens = apply(dens2D,1,sum)

    # Renormalize
    dens = normalize_density(dens,mids)

    out = list(x=mids,y=dens)
    return(out)
}

#' @export
combine <- function(x,y,normalize=FALSE) 
{ # Function to combine all combinations of two vectors
  
  tmp <- expand.grid(x,y)
  out <- tmp[[1]]*tmp[[2]]
  
  if (normalize) out <- out/sum(out)
  
  cat("Combined values:",length(out),"\n")
  
  return(out)
}

#' @export
plot.histweights <- function(dist)
{
  col <- rep(2,length(dist$x))
  col[dist$i95] <- "violet"
  col[dist$i66] <- 1
  col[dist$i10] <- 4
  plot(dist$x,dist$y,col=col,xlab="Temp (°C)",ylab="Weighting",xlim=c(0,6))

}

#' @export
pointfit2 <- function(x0,y0,col=1,lwd=1,pch=1,lty=1)
{
  fit <- lm(y0~x0+I(x0^2))
  x1 <- sort(x0); y1 <- predict(fit,newdata=data.frame(x0=x1))
  points(x0,y0,pch=pch,col=col)
  lines(x1,y1,col=alpha(col,60),lwd=lwd,lty=lty)

}
