 
## This file contains general functions needed
## for statistical analysis of individual
## time series


### NON-LINEAR TRENDS ####
require(Rssa)

make_gaps <- function(y,per=0.2)
{
  n = length(y)
  y[sample(1:n,size=round(n*0.1))] = NA 
  y[1:round(n*per)] = NA 
  y[1:round(n*per)+round(n*0.4)] = NA
  y[round(n*0.9):(n)] = NA
  return(y)
}

get_gaps <- function(y,gapmax=10,consecmin=5)
{
  
  # First populate temporary vector/array
  y1 = y 

  # If it's a 2D array (nmonths x ntime), collapse to 1 dim
  if (!is.null(dim(y1))) y1 = apply(y1,2,mean)
  
  # Get length of original (collapsed) series
  ny = length(y1)

  # Add one to each side to make sure ends are always 'good'
  y1 = c(0,y1,0)
  
  # Get length of padded series
  n = length(y1)
  
  # Determine valid points and the differences between indices
  good  = which(!is.na(y1))
  dgood = c(diff(good),1)
  diffs = rep(0,n)
  diffs[good] = dgood
  
  # Breaks occur when the gap in indices is bigger than the
  # maximum
  breaks = which(abs(diffs) > gapmax)
  nbreaks = length(breaks)

  # If breaks exist, separate time series into segments
  if ( nbreaks > 0 ) {
    
    # Store the indices of each set of breaks
    bb = array(NA,dim=c(2,nbreaks))
    for (q in 1:nbreaks) {
      b0 = breaks[q]
      tmp = good
      tmp[tmp<=b0] = NA
      b1 = tmp[ which(!is.na(tmp))[1] ]
      
      # Make sure to adjust for extra boundary points
      # to get the right indices of actual vector
      if (b0==1) b0 = b0+1
      if (b1==n) b1 = b1-1
      bb[,q] = c(b0,b1)-1
    }
    
    # Now add boundary breaks
    bb = cbind(c(0,0),bb,c(n,n))
    
    # Get indices of actual segments of interest
    # Store each segment by number
    x = c(1:length(y1))
    segments = x*0
    k = 0
    for (q in 2:dim(bb)[2]) {
      tmp = which(x %in% c(bb[2,q-1]:bb[1,q]))
      if (length(tmp)>consecmin) {
        k = k+1
        segments[tmp] = k
      }
    }
    
    # Limit segments back to original series
    segments = segments[2:(ny+1)]

  } else {
    segments = rep(1,ny)
  }
  
  return(segments)
}

test_gaps <- function(y)
{
  y0 = y
  y = make_gaps(y0)
  x = c(1:length(y))
  
  seq1 = get_gaps(y)

  plot(x,y0,type="n")
  abline(v=x[which(seq1>0)],col=alpha(2,50))
  points(x,y0,col='grey70')
  points(x,y)
}

# Just return y values from approx
approxy <- function(...) approx(...)$y

series_fill <- function(y)
{ # Fill in the missing values of the time series 
  # with a linear interpolation
  
  twoD = FALSE
  if (!is.null(dim(y))) twoD = TRUE 

  # Generate a generic x vector
  # and a new y vector
  x = c(1:length(y))
  if (twoD) x = c(1:dim(y)[2])
  ynew = y
  
  # Find indices of missing values (1D representation)
  ii = which(is.na(y))

  # Approximate the missing values with a linear fit
  if (length(ii)>0) {
    if (!twoD) {
      ynew[ii] = approxy(x,y,xout=x[ii],rule=2)
    } else {
      ynew = t( apply(y,1,approxy,xout=x,rule=2) )
    }
  }

  return(list(y=ynew,ii=ii))
}

max_period <- function(y)
{ # Returns the period of maximum spectral power of a time series (or several)
  
  sp = spectrum(y,plot=F)
  k = which.max(sp$spec)

  pmax = 1/sp$freq[k]
  return(pmax)
}

add_buffer.internal <- function(y,nrep=1,nbuff=floor(length(y)*0.2),scale=0.5,type="mean")
{
  ### Modify the time series to facilitate spectral analysis over the whole time series:
  # Preferred Method: add mirrored and inverted nbuff points to each end of y, remove them later
  ny = length(y)
  
  # Get a window of the beginning and end of our time series
  y00a = y[1:nbuff]
  y11a = y[(ny-nbuff+1):ny]
  
  if (type=="mean") {
    # Make buffer points equal mean plus random noise
    y00 = rep(mean(y00a),nbuff*nrep) + rnorm(nbuff*nrep,sd=sd(y00a)*scale)
    y11 = rep(mean(y11a),nbuff*nrep) + rnorm(nbuff*nrep,sd=sd(y11a)*scale)
    y1  = c(y00,y,y11)
  
  } else if (type=="linear") {
    # Make buffer points equal to linear trend plus random noise
    x   = c(1:nbuff)
    x0  = c(-(nbuff*nrep-1):0)
    x1  = c(1:(nbuff*nrep))+nbuff
    y00 = as.vector(predict(lm(y00a~x),newdata=list(x=x0))) + rnorm(nbuff*nrep,sd=sd(y00a)*scale)
    y11 = as.vector(predict(lm(y11a~x),newdata=list(x=x1))) + rnorm(nbuff*nrep,sd=sd(y11a)*scale)
    y1  = c(y00,y,y11)
  }

  # # Get initial and final index of actual time series for later use
  # i0 <- nrep*nbuff + 1
  # i1 <- ny + i0 - 1

  #return(list(y=as.vector(y1),i0=i0,i1=i1)) #,y00a=y00a,y11a=y11a,x0=x0,x1=x1))

  return(y1)
}

add_buffer <- function(y,nrep=1,nbuff=floor(length(y)*0.2),scale=0.5,type="mean")
{ # Now can handle 2D time series (nm X nt)
  ### Modify the time series to facilitate spectral analysis over the whole time series:
  # Preferred Method: add mirrored and inverted nbuff points to each end of y, remove them later
  y0 = y
  if (is.null(dim(y))) dim(y0) = c(1,length(y0))
  ny = dim(y0)[2]

  y1 = t( apply(y0,1,add_buffer.internal,nrep=nrep,nbuff=nbuff,scale=scale,type=type) )
  
  # Get initial and final index of actual time series for later use
  i0 <- nrep*nbuff + 1
  i1 <- ny + i0 - 1
  
  # Make sure output dimensions are equivalent
  if (is.null(dim(y))) y1 = as.vector(y1)

  return(list(y=y1,i0=i0,i1=i1)) #,y00a=y00a,y11a=y11a,x0=x0,x1=x1))
}

test_buffer <- function(y,nbuff=13,nrep=2)
{
  ny = length(y)
  test = list(x=c(1:length(y)),y=y)
  ybuf = add_buffer(y,nbuff=nbuff,nrep=nrep,type="linear")
  i0 = ybuf$i0 - nrep*nbuff
  i1 = ybuf$i1 - nrep*nbuff
  x0 = ybuf$x0
  x1 = ybuf$x1
  y00a = ybuf$y00a
  y11a = ybuf$y11a

  x  = c(1:nbuff)
  xtmp = c((1-nbuff*nrep):(ny+nbuff*nrep))
  plot(xtmp,ybuf$y,col="grey50")
  points(test)
  abline(v=c(i0,i1))
  lines(c(x0,x),predict(lm(y00a~x),newdata=list(x=c(x0,x))),col=2)
  lines(c(x,x1)+(ny-nbuff),predict(lm(y11a~x),newdata=list(x=c(x,x1))),col=2)
}

ssatrend.auto <- function(y,pmin=30,frequency=12,pbuff=0,L=NULL,fill=TRUE,smooth=TRUE,
                          Fmax=8,old=F)
{ # Compute the trend from SSA decomposition
  # Chooses close to maximum L (a little less than 1/2 the length of y)
  # Then limits the trend to the eigenvectors above a certain period of interest
  # (ie the long-term trend)
  
  # Save dimensions of y and vectorize (in case of monthly values)
  dims = dim(y)
  y1   = as.vector(y)
  ny   = length(y1)
  i0   = 1
  i1   = ny 
  
  # Get dy (time step)
  dy = 1/frequency 

  # Get L (rep of pmin, but less than ny/2, and greater than zero!)
  if (is.null(L)) L = (ny/2) - ( (ny/2) %% floor( min(pmin/dy,ny/2)) )

  # How many repetitions of buffer points?
  nbuff = min( c(floor(ny/2),floor(pbuff/dy)) )
  nrep  = 1

  # For a very short time series, make sure
  # buffer points will be added for stability
  if (ny*dy < pmin) {
    nbuff = ny
    nrep  = ceiling( pmin/(ny*dy) )
  }

  # If we are using old method, set key parameters here
  if (old) {
    L      = floor(15/dy)  # 15 year window
    Fmax   = 1             # Only the first eigenvector
    nrep   = 0             # No buffer points
    smooth = FALSE         # Don't use preliminary ssa smoothing
  }

  ## If fill requested, make sure that NAs are filled in
  if (fill) { 
    tmp   = series_fill(y1)
    ii.na = tmp$ii 
    y1    = tmp$y 
  }
  
  ## Add buffer points if desired
  if (nbuff > 0) {
    tmp = add_buffer(y1,nrep=nrep,nbuff=nbuff)
    y1  = tmp$y 
    i0  = tmp$i0
    i1  = tmp$i1
  }
  
  # First smooth the trend via a very small window
  # to elimate any noise problems
  if (smooth) {
    ssa0 = ssa(y1,L=1/dy)
    y1 = reconstruct(ssa0,groups=list(c(1)))$F1
  }

  # Get a new ssa object and reconstruct the first 12 eigen vectors
  ssa = ssa(y1,L=L)
  new = reconstruct(ssa,groups=as.list(c(1:Fmax)))
  
  # Convert list of eigenvectors to a data.frame of time series
  new = as.data.frame(new)
  
  # Determine which eigenvectors show periods greater than pmin
  qq = which( apply(new,2,max_period)*dy >= pmin )

  # Sum the eigen vectors of interest to get the trend
  trend = apply(as.data.frame(new[qq]),1,sum)
  
  ### Remove the extra boundary buffer points from the series
  trend = trend[i0:i1]
  
  # Now put missing values back in 
  if (fill) trend[ii.na] = NA 
  
  # If y is a time series make sure output is a time series!  
  if (is.ts(y)) trend = ts(trend,start=start(y),end=end(y),frequency=frequency(y))

  # Make sure dims are the same as initial data 
  dim(trend) = dims 
  
  # That's it! Return the trend...
  return(trend)
}

ssatrend.fill <- function(y,pmin=30,frequency=12,L=NULL,pbuff=0,smooth=T,
                          gapmax=(pmin*frequency)*0.5,consecmin=(10*frequency))
{ # Perform ssatrend by consecutive segments of real data in a time series
  
  trend = y*NA 

  if ( sum(is.na(y)) < length(y) ) {
    # Determine the segments of the series
    segments = get_gaps(y,gapmax=gapmax,consecmin=consecmin)
    
    ss = unique(segments[segments>0])
    
    # Loop over segments and get trend for each segment
    for (s in ss) {
      ii = which(segments == s)
      y1 = y[ii]
      trend1 = ssatrend.auto(y1,pmin=pmin,frequency=frequency,pbuff=pbuff,L=L,fill=TRUE,smooth=smooth)
      trend[ii] = trend1
    }
  
  }
  
  if (is.ts(y)) trend = ts(trend,start=start(y),end=end(y),frequency=frequency(y))

  return(trend)
}

ssatrend.simple <- function(y,L=15,eigen=1,fill=TRUE,nbuff=11,buffer="linear",svd_method=NULL)
{ # Compute the trend from SSA decomposition
  
  # Save dimensions of y and vectorize (in case of monthly values)
  dims = dim(y)
  y1   = as.vector(y)
  ny   = length(y1)
  i0   = 1
  i1   = ny 

  # For a very short time series, make sure
  # buffer points will be added for stability
  nrep  = 1
  if (ny < L) {
    nbuff = ny
    nrep  = 2*ceiling( L/ny )
  }

  # If fill requested, make sure that NAs are filled in
  if (fill) { 
    tmp   = series_fill(y1)
    y1    = tmp$y 
    ii.na = tmp$ii 
    if (sum(is.na(y1))>0) cat("Error: NA present!",y1,"\n")
  }
  
  ## Add buffer points if desired
  if (nbuff > 0) {
    tmp = add_buffer(y1,nrep=nrep,nbuff=nbuff,type=buffer)
    y1  = as.vector(tmp$y)
    i0  = tmp$i0
    i1  = tmp$i1
  }
  
  # Get a new ssa object and
  # Reconstruct the first several eigen vectors
  ssa = ssa(y1,L=L,svd_method=svd_method)
  new = reconstruct(ssa,groups=as.list(eigen))
  
  # Sum the eigen vectors of interest
  trend = apply(as.data.frame(new[eigen]),1,sum)
  
  ### Remove the extra boundary buffer points from the series
  trend = trend[i0:i1]

  # Now put missing values back in 
  if (fill) trend[ii.na] = NA 
  
  # Make sure dims are the same as initial data 
  dim(trend) = dims 
  
  if (is.ts(y)) trend = ts(trend,start=start(y),end=end(y),frequency=frequency(y))

  # That's it! Return the trend...
  return(trend)
}

ssatrend2D.simple <- function(y,L=15,eigen=1,fill=TRUE,nbuff=11,buffer="linear",svd_method=NULL)
{ # Compute the trend from SSA decomposition
  
  # Save dimensions of y and vectorize (in case of monthly values)
  dims = dim(y)
  y1   = y
  ny   = dims[2]
  i0   = 1
  i1   = ny

  # For a very short time series, make sure
  # buffer points will be added for stability
  nrep  = 1
  if (ny < L) {
    nbuff = ny
    nrep  = 2*ceiling( L/ny )
  }

  # If fill requested, make sure that NAs are filled in
  if (fill) { 
    tmp   = series_fill(y1)
    y1    = tmp$y 
    ii.na = tmp$ii 
    if (sum(is.na(y1))>0) cat("Error: NA present!",y1,"\n")
  }
  
  ## Add buffer points if desired
  if (nbuff > 0) {
    tmp = add_buffer(y1,nrep=nrep,nbuff=nbuff,type=buffer)
    y1  = tmp$y
    i0  = tmp$i0
    i1  = tmp$i1
  }
  
  # Get a new ssa object and
  # Reconstruct the first several eigen vectors
  # svd_method = NULL, let ssa decide
  # svd_method = "propack", for small values of L 
  ssa = ssa(y1,kind="2d-ssa",L=c(1,L),neig=2,svd_method=svd_method)
  new = reconstruct(ssa,groups=as.list(1,eigen))
  
  # Sum the eigen vectors of interest
  trend = new$F1
  #trend = apply(as.data.frame(new[eigen]),1,sum)
  
  ### Remove the extra boundary buffer points from the series
  trend = trend[,i0:i1]

  # Now put missing values back in 
  if (fill) trend[ii.na] = NA 
  
  # Make sure dims are the same as initial data 
  dim(trend) = dims 
  
  if (is.ts(y)) trend = ts(trend,start=start(y),end=end(y),frequency=frequency(y))

  # That's it! Return the trend...
  return(trend)
}

ssatrend.fill.simple <- function(y,L=15,nbuff=11,buffer="linear",svd_method=NULL,
                                 gapmax=L,consecmin=L+5,twoD=FALSE)
{ # Perform ssatrend by consecutive segments of real data in a time series
  
  trend = y*NA 

  # If there are some data points available, peform ssatrend
  if ( sum(is.na(y)) < length(y) ) {
    
    # Determine the segments of the series
    segments = get_gaps(y,gapmax=gapmax,consecmin=consecmin)
    
    ss = unique(segments[segments>0])
    
    if (!twoD) {
      # Loop over segments and get trend for each segment
      for (s in ss) {
        kk = which(segments == s)
        y1 = y[kk]
        trend1 = ssatrend.simple(y1,L=L,nbuff=nbuff,buffer=buffer,svd_method=svd_method)
        trend[kk] = trend1
      }
    } else {
      # Loop over segments and get trend for each segment
      for (s in ss) {
        kk = which(segments == s)
        y1 = y[,kk]
        trend1 = ssatrend2D.simple(y1,L=L,nbuff=nbuff,buffer=buffer,svd_method=svd_method)
        trend[,kk] = trend1
      }
    }
  }
  
  if (is.ts(y)) trend = ts(trend,start=start(y),end=end(y),frequency=frequency(y))

  return(trend)
}

# Use the wrapper ssatrend.fill as the default method!
ssatrend = ssatrend.fill.simple

### END NON-LINEAR TRENDS ####

## Spectra
my.spectrum <- function(y,...,x.labs=c(5,10,15,20,30,50,100,500))
{
  sp = spectrum(y,plot=F)
  plot(sp$freq,sp$spec,type='n',ann=F,axes=F,xlim=range(1/x.labs),...)
  mtext(side=1,line=2,"Period (a)")
  mtext(side=2,line=2,las=0,"Spectral power")
  x.at = 1/x.labs
  axis(1,at=x.at,labels=x.labs)
  axis(2)

  lines(sp$freq,sp$spec)
  
  return(sp)
}

## DIM'S C++ FUNCTIONS (translated)

# Error function (and complementary error function)
erf  <- function(x) 2 * pnorm(x*sqrt(2)) - 1
erfc <- function(x) 2 * pnorm(x*sqrt(2), lower = FALSE)

pdf_dim <- function(x,t,sigma,mu1,mu0 )
{ # Intermediate function for p_x()
  pdf0 = 1/sqrt(2*pi*sigma^2) * exp( -(x-mu0-mu1*t)^2 / (2*sigma^2) )

  return(pdf0)
}

cdf <- function(x,t,sigma,mu1,mu0 )
{ # Intermediate function for p_rec()
  f_x = 0.5 - 0.5 * erf( (x-mu0-mu1*t)/(sigma*sqrt(2)) )
  return(f_x)
}

p_x <- function(x,t,sigma,mu1,mu0,dx=1e-8) 
{ # Intermediate function for p_rec() 
  if ( t < 1 ) stop("t must be >= 1.")
  ff = 1
  if (t>1) {
    for ( ti in 1:(t-1) ) ff = ff * (1.0 - cdf(x,ti,sigma,mu1,mu0))
  }

  # p = (cdf(x,t,sigma,mu1,mu0)-cdf(x+dx,t,sigma,mu1,mu0))*ff / dx
  p = pdf_dim(x,t,sigma,mu1,mu0)*ff

  return(p)
}

p_rec <- function(t,sigma=1,mu1=0,mu0=0,model=0,subdivisions=1e2)
{ # Calculate the probability of a record at timestep n
  # model==0: 1/n
  # model==1: Integration
  # model==2: Linear approximation
  #            **Only valid for large n, so to avoid errors, n < 6 will not be permitted

  if ( model==2 ) {
        
    p_rec = alpha = t*NA

    #if ( min(t) < 6 ) warning("t < 6 will not be calculated (set to NA)")

    ii <- which(t > 6)
    alpha[ii] = 2*sqrt(pi)/exp(2) * sqrt(log(t[ii]^2/(8*pi)))

    p_rec = 1/t + alpha * mu1/sigma

    
  } else if ( model==1 ) {
    
    # In case calculating for more than 1 value of t, use a loop:
    p_rec = t*NA
    
    for (i in 1:length(p_rec) ) {
      ## For integration, use the 1-d integration function:
      xrange = c(-20*sigma,20*sigma)+mu1*t[i]
      #xrange = c(-Inf,Inf)
      tmp = integrate(p_x,xrange[1],xrange[2],subdivisions=subdivisions,t=t[i],sigma=sigma,mu1=mu1,mu0=mu0)
      p_rec[i] = tmp$value
    }
  
  } else if ( model==0 ) {

    p_rec = 1/t

  }
  
  # Set values below zero to a very small value
  # (not zero that divisions are ok..)
  p_rec[p_rec<=0] = 1e-8

  return(p_rec)
}

gen.lookup <- function(t=c(1:200),mu1sigma=seq(-0.2,0.2,by=0.01))
{ # Generate a lookup table for t values of p_rec and various mu1/sigma ratios
  
  if (file.exists("lu.Rdata") ) {
    load("lu.Rdata")
  } else {
    lu = list(mu1sigma=round(mu1sigma,2))
    lu$prec1 = lu$prec2 = list()
    for (i in 1:length(lu$mu1sigma)) {
      lu$prec1[[i]] = p_rec(t,sigma=1,mu1=lu$mu1sigma[i],model=1)
      lu$prec2[[i]] = p_rec(t,sigma=1,mu1=lu$mu1sigma[i],model=2)
    }
    save(lu,file="lu.Rdata")
  }

  return(lu)
}

p_rec.lookup <- function(t,sigma=1,mu1=0,mu0=0,model=1)#,lu)
{ # Extract p_rec series of interest from a pre-calculated look-up table,
  # Find the time series of expected records that matches the desired mu1sigma ratio
  
  p_rec = t*NA
  if ( !is.na(sigma+mu1+mu0) ) {
    i = which(lu$mu1sigma==round(mu1/sigma,2))
    #cat("mu1sigma=",round(mu1/sigma,2),"\n")
    if ( model==1 ) p_rec = lu$prec1[[i]][t]
    if ( model==2 ) p_rec = lu$prec2[[i]][t]
  }

  return(p_rec)
}

Rx <- function(t,ntot=10,sigma=1,mu1=0,mu0=0,model=0) 
{ # Calculate the expected number of records for a given time interval
  m = ntot-1
  # Rx = 0
  # for ( ti in (t-m):t ) Rx = Rx + p_rec(ti,sigma=sigma,mu1=mu1,mu0=mu0,model=model)
  
  Rx = sum( p_rec(c((t-m):t),sigma=sigma,mu1=mu1,mu0=mu0,model=model) )

  return(Rx)
}

Xx <- function(t,ntot=10,sigma=1,mu1=0,mu0=0,model=1) 
{ # Determine the ratio of records, calculating them at each time first
  X1 = Rx(t,ntot,sigma,mu1,mu0,model=model)
  X0 = Rx(t,ntot,sigma,mu1,mu0,model=0)

  return(X1/X0)
}

Xrec <- function(p1,p0,ii)
{ # Determine the ratio of records over a given time interval (indices ii)
  X1 = sum(p1[ii])
  X0 = sum(p0[ii])

  return(X1/X0)
}

## NORMALITY TESTS ####
## (with simplified output for analysis)

my.ks.test <- function(x,y,...)
{ # Perform ks.test, only return D-statistic.
  tmp = ks.test(x,y,...)
  D = as.numeric(tmp$statistic)
  p = as.numeric(tmp$p.value)
  return(D)
}

my.bptest <- function(...)
{ # Perform bp.test, only return BP statistic.
  #require(lmtest)
  tmp = bptest(...)
  BP  = as.numeric(tmp$statistic)
  p   = as.numeric(tmp$p.value)
  return(BP)
}

## END NORMALITY TESTS ####

## NOW TIME SERIES FUNCTIONS ####

series_gen <- function(n=100,a=0.1,missing=integer(0))
{ # Generate a random series with linear trend
  # for testing

  x = c(1:n)
  y = rnorm(length(x)) + x*a
  y[missing] = NA
  trend = ssatrend(y,L=5)
  sd = sd(y-trend)
  
  return(list(x=x,y=y,trend=trend,sd=sd))
}

series_extremes <- function(y)
{ # Calculate when extremes occur in time series y
  
  # Generate empty arrays to hold extremes for each time step
  exhi = y*0
  exlo = y*0

  nt = length(y)

  # Temperature level reached by each realization (start at -1000 degC to ensure first value is extreme)
  exhival = -1e3
  exloval =  1e3
  
  # Loop over each time step and determine if each realization has an extreme
  # (as well as for the original data, and the detrended data)
  for ( k in 1:nt ) {
    
    # Update the extremes vector to reflect the new found levels
    # Add to new extremes to counts
    
    if (!is.na(y[k])) {

      # Upper extremes
      if ( y[k] > exhival ) {
        exhival    = y[k]
        exhi[k]    = 1
      }
      
      # Lower extremes
      if ( y[k] < exloval ) {
        exloval    = y[k]
        exlo[k]    = 1
      }
    }
    
  }

  return(data.frame(exhi=exhi,exlo=exlo))
}

series2D_extremes <- function(yy)
{ # Calculate when extremes occur in 2D time series yy (nm X n)
  # (To get low extremes, input negative time series (-yy) )

  # Get dimensions of input
  nm = dim(yy)[1]
  nt = dim(yy)[2]

  # Generate empty arrays to hold extremes for each time step
  exhi = array(0,dim=dim(yy))

  # Temperature level reached by each realization (start at -1000 degC to ensure first value is extreme)
  exhival = rep(-1e3,nm)
  
  # Loop over each time step and determine if each realization has an extreme
  # (as well as for the original data, and the detrended data)
  for ( k in 1:nt ) {
    
    # Update the extremes vector to reflect the new found levels
    # Add new extremes to counts
    
    # Upper extremes
    ii = which(yy[,k] > exhival)
    if (length(ii)>0) {
      exhival[ii] = yy[ii,k]
      exhi[ii,k]  = 1
    }
    
  }

  return(exhi)
}

series_stats <- function(x,y,fit="ssa",pmin=30,L=NULL,nmin=5,p_rec_func=p_rec.lookup)
{
  nt    = length(y)
  x0    = x-x[1]+1

  # How many points in the series are good?
  ngood = length(which( !is.na(y) ))
  
  ## Define ALL output values
  trend = y0 = prec0 = prec1 = prec2 = y*NA
  sigma = mu1 = mu0 = ksD = BP = NA
  sigmas = rep(NA,5)
  ex = data.frame(exhi=y*NA,exlo=y*NA)

  ## Only perform analysis if there are enough data points
  if ( ngood >= nmin ) {

    # 1. Get trend and trend parameters
    # if ( fit == "linear" ) {
    #   tmp   = lm(y~x0)
    #   trend = as.vector(predict(tmp,newdata=data.frame(x0=x0)))
    #   mu0   = as.numeric(tmp$coefficients[1])
    #   mu1   = as.numeric(tmp$coefficients[2])

    # } else {

      trend = ssatrend.fill(y,L=15)
    # }
    
    # 2. Get sigma of detrended data
    y0    = y - trend
    sigma = sd(y0,na.rm=T)

    # 3. Check Gaussinity of detrended data
    # (if statistic "D" is small, likely Gaussian)
    # (also, if p.value is big, more confidence)
    # if ( !is.na(sum(y0)) ) {
    #   # tmp = ks.test(x=y0,y="pnorm")
    #   # y0gauss = c(tmp$statistic,tmp$p.value)
    #   ksD = my.ks.test(x=y0,y="pnorm")
    #   #if ( fit == "linear" ) BP  = my.bptest(y~x0)
    # }

    # # If linear fit was used, check how many theoretical distribution of extremes
    # prec0 = 1/x0  #p_rec(t=c(1:nt),model=0)
    # if ( fit == "linear" ) {
    #   prec1 = p_rec_func(t=x0,sigma=sigma,mu1=mu1,mu0=mu0,model=1)
    #   prec2 = p_rec_func(t=x0,sigma=sigma,mu1=mu1,mu0=mu0,model=2)  
    # }

    # 3. Get extremes
    ex = series_extremes(y)
  }
  
  return(list(x=x,y=y,trend=trend,y0=y0,sigma=sigma,exhi=ex$exhi,exlo=ex$exlo))
}

get.intervals <- function(x,weights=NULL,norm=FALSE)
{ # Function to determine confidence/credence intervals
  # from random samples

  # Number of points
  n = length(x)
  
  if (is.null(weights)) weights = rep((1/n),n) 

  # Order the samples from smallest to largest
  i1 = order(x)
  x1 = x[i1]
  weights1 = weights[i1]
  
  # Determine the intervals that we will output
  # Should be c(2sigma-,1sigma-,mean,1sigma+,2sigma+)
  conf.ranges = c(2.3,15.9,50.0,84.1,97.7)

  # Cumulatively sum the weights to get confidence intervals
  intervals = cumsum(weights1*100)
  
  # Generate an output vector with the same length as 
  # conf.ranges
  out = conf.ranges*NA

  # Determine the limits of each confidence interval
  # and store in output vector
  for ( q in 1:length(out) ) {
    ii = which( intervals <= conf.ranges[q] )
    out[q] = max(x1[ii])
    
    if (norm) out[q] = sum(x1[ii]*weights1[ii])
  }

  # Return the limit of each confidence interval
  return(out)
}

get.intervals.binom <- function(x)
{ # Function to determine confidence/credence intervals
  # from random samples from a binomial distribution
  
  success = sum(x)
  total   = length(x)

  conf.levels = c(0.682,0.954)
  
  out = rep(NA,5)

  tmp = binom.test(success,total,conf.level=conf.levels[1])
  out[c(2,4)] = as.numeric(tmp$conf.int)
  tmp = binom.test(success,total,conf.level=conf.levels[2])
  out[c(1,5)] = as.numeric(tmp$conf.int)

  out[3] = as.numeric(tmp$estimate)
  
  # Return confidence intervals
  # Should be c(2sigma-,1sigma-,mean,1sigma+,2sigma+)
  return(out)
}


series_mc <- function(trend,sigma,n=100)
{ # Function to perform monte carlo experiment
  # given a trend (vector the length of time series), and
  # sigma - the standard deviation of the time series 
  # n is the number of Monte Carlo samples desired
  
  # Determine the length of the time series
  nt = length(trend)

  # If the stand. dev. is a vector, then
  # replicate it to be the size of the output
  # Monte Carlo array
  if (length(sigma)>1) sigma = rep(sigma,n)

  # Generate all random realizations (length of time series by number of realizations)
  mc = array(rnorm(nt*n,sd=sigma,mean=0),dim=c(nt,n))
  
  # Add the trend onto each realization (as an array for efficiency)
  mct = mc + array(rep(trend,n),dim=c(nt,n))
  
  # Temperature level reached by each realization (start at -1000 degC to ensure first value is extreme)
  exhival = numeric(n) - 1000.0
  exloval = numeric(n) + 1000.0
  
  # Extremes for each time step
  exhi = mct*0
  exlo = mct*0
  
  # Confidence interval output arrays
  # c(2sigma-,1sigma-,mean,1sigma+,2sigma+)
  confhi = conflo = as.data.frame(array(NA,dim=c(nt,5)))
  names(confhi) = names(conflo) = c("sig2a","sig1a","mid","sig1b","sig2b")

  # Loop over each time step and determine if each realization has an extreme
  # (as well as for the original data, and the detrended data)
  for ( k in 1:nt ) {
    
    # Update the extremes vector to reflect the new found levels
    # Add to new extremes to counts
    # Calculate the confidence ranges 

    # Upper extremes
    ii1  = which( mct[k,] > exhival )
    exhival[ii1] = mct[k,ii1]
    exhi[k,ii1]  = 1
    confhi[k,]   = get.intervals.binom(exhi[k,])

    # Lower extremes
    ii2 = which( mct[k,] < exloval )
    exloval[ii2] = mct[k,ii2]
    exlo[k,ii2]  = 1
    conflo[k,]   = get.intervals.binom(exlo[k,])
  }
  
  # Return the estimated conf ranges for exhi and exlo
  return(list(exhi.mc=confhi,exlo.mc=conflo))   # exhi=exhi,exlo=exlo,mct=mct
}


## series_seasmean 
series_seasmean <- function(var12,months=NULL)
{ # Get the seasonal average, making sure to shift december to next year
  nt = dim(var12)[2]
  var12[12,] = c( NA, var12[12,1:(nt-1)] )

  
  s = list()
  if (!is.null(months)) {
    s[[1]] = months
  } else {
    s[[1]] = c(12,1,2)
    s[[2]] = c(3,4,5)
    s[[3]] = c(6,7,8)
    s[[4]] = c(9,10,11)
  }
  
  ns = length(s)

  vars = array(NA,dim=c(ns,nt))

  for (q in 1:ns) {
    vars[q,] = apply(var12[s[[q]],],MARGIN=2,FUN=mean)
  }

  return(vars)
}

