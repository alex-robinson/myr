## FILTERS ##



ssatrend.jorge <- function(x,window=length(x)/2,Tmin=0,Tmax=1e6) 
{ # Obtain the singular-sprectrum trend of a time series
  # (not finished yet!)

  x.decomp = decompSSA(x,window)
  
  freqs = x.decomp$freq
  eigen = x.decomp$U
  
  # Output the list of frequencies available
  # To know the eigenvalue of the lowest frequency
  i <- which.min(freqs)

# Now we search the eigenvectors for the one representing the trend.
  e = x.decomp$U
  indx = which.max(abs(sum(e))/sum(abs(e)))
# Just need the one vector representing the trend
  e = e[,indx]
# Now create a filter using this vector
  efilt = convolve(e,flipud(e),type="o")
# Normalize
  efilt = efilt/sum(efilt)
  
# Now we have to pad the original series so that we can get a full trend.
# We do this by adding values to the start and end based on the linear trend
# of the first window of values.
  n = length(x)
  df = data.frame(x=x[1:window],time=seq(1,window))
  res = lm(x~time,data=df)
  left.df = data.frame(time=seq(1-window,0))
  x.pad = c(predict.lm(res,left.df),x)
  df = data.frame(x=x[(n-window+1):n],time=seq(1,window))
  res = lm(x~time,data=df)
  right.df = data.frame(time=seq(window+1,2*window))
  x.pad = c(x.pad,predict.lm(res,right.df))
  
# Now just run it through the filter
  result = filter(x.pad,efilt)
  result = result[(window+1):(length(result)-window)]
  
  return(result)
}


#Hamming and Blackman windows from 0..M
Hamming   <- function( M, x=c(0:M)/M ) { 0.54-0.46*cos(2*pi*x) }
Blackman  <- function( M, x=c(0:M)/M ) { 0.42-0.5*cos(2*pi*x)+0.08*cos(4*pi*x) }

#sinc function of frequency f
sinc      <- function( M, fc, x=c(0:M)-M/2 ) { ifelse( x==0, 2*pi*fc, sin(2*pi*fc*x)/x ) }

wlpf_internal <- function(y, dt, fc, Mmax )
{ # Internal function to calculate lopass filter,
  # Should only be called from wrapper wlpf

  # Length of time series
  ny    <- length(y)
  
  # Determine the width of the sinc kernel
  # (Limited by Mmax if desired)
  M  <- floor(ny/2)
  if ( !is.na(Mmax) ) M <- min(M, Mmax)
  
  # Bandwidth
  BW <- 4/M
  
  ### Output filter information to screen
  cat("fc=",fc,"; pc=",1/fc,"; dt=",dt,"; Filter bandwidth=",BW,
      "; Window half-width M=",M,"\n")

  ### Modify the time series to facilitate spectral analysis over the whole time series:
  # Preferred Method: add mirrored m+1 points to each end of y, remove them later
  y00 <- y[1:(M/2+1)]
  y11 <- y[(ny-M/2):ny]  #*0 + mean(y)
  y1  <- c(rev(y00),y,rev(y11))
  
  # Get initial and final index of interest for later use
  i0 <- (M/2+1)+(M/2+1)+1
  i1 <- ny + i0 - 1
  
  # Other Methods: 
  # (1) leave time series the same, lose (M/2+1) points at the end of series
  # (2) Only add new points onto the end (may cause strange behavior at beginning of series)
  #y1   <- y
  #y1  <- c(y,rev(y11))
  
  # Get mean and sd, normalize  data
  ymean <- mean(y)
  ysd   <- sd(y)
  y1    <- (y1-ymean)/ysd
  
  ### Generate the sinc kernel and Blackman window. Normalize.
  rk1 <- sinc( M, fc*dt ) 
  bk1 <- rk1 * Blackman(M)
  bk1 <- bk1/sum(bk1)
  bk  <- bk1
  
  ### Convolve time series with the filter kernel
  # using FAST FOURIER TRANSFORM

  # Pad the filter with zeros
  k   <- c(bk, rep(0,length(y1)-length(bk)))

  # Convolve data and filter kernel
  fy  <- fft(fft(k)*fft(y1), inverse=TRUE)
  out <- Re(fy)

  # Rescale to match original data
  out <- (out-mean(out))/sd(out)
  out <- out*ysd + ymean
  
  ### Remove the extra points from the series
  ### (to account for the extra buffer points we added,
  ### plus the initial half-width to be thrown away)
  out <- out[i0:i1]

  return(out)
}

freqfilter <- function( x=NA, y, dt=NA, pc, fc=1/pc, Mmax=NA, bp=FALSE, vec=TRUE )
{ # A wrapper for calls to the window-sinc filter 
  ## Windowed sinc low pass filter
  #   y  - vector to filter
  #   dt - time interval between measurements (s)
  #   fc - low pass frequency (Hz)
  
### Organize the input data, make sure it is consistent
  # If x exists,determine dt (assumes equal spacing of x)
  if (!is.na(x[1])) dt = x[2]-x[1]
  # If dt still hasn't been created, produce an error
  if (is.na(dt)) {
    stop("Please either specify the x vector or a dt value.\n")
  }
  
  if (bp & length(fc) != 2) {
    stop("For a bandpass filter, please specify a range of frequencies/periods.\n")
  }
  
  # For bandpass filter, make sure frequencies are calculated in ascending order
  if (bp) fc = sort(fc)
  
  out = array(NA,dim=c(length(fc),length(y)))

  for (q in 1:length(fc)) out[q,] = wlpf_internal(y=y,dt=dt,fc=fc[q],Mmax=Mmax)
  if (length(fc)==1) out = as.vector(out)

  # For bandpass filter, subtract lo-pass filter from high-pass filter
  if (bp) out.bp = out[2,] - out[1,]
  
  if (!vec) {
    out = list(pc=pc,fc=fc,yf=out,bp=out.bp)
  }

  if (vec & bp) out = out.bp 
  
  # Output the frequency-filtered vector
  return(out)
}

