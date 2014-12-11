

## LOADING BINARY (GRADS) FILES


## Load 1d binary file
#' @export
load_grads1D <- function(fldr="/scratch/01/andrey/PPSCM",ctl="grads/history.ctl",file="OUT/history.dat",n=NULL)
{
  ## Load variable names from CTL file
  ctl <- scan(file.path(fldr,ctl),what="character",sep="\n")
  ctl <- strsplit(ctl,split=" ")
  nms <- vector(mode="character")
  vars <- FALSE
  nt = NULL 
  for ( i in 1:length(ctl)) {
    s <- ctl[[i]][1]
    if ( s == "ENDVARS" )vars <- FALSE
    if ( vars ) nms <- c(nms,s)
    if ( s == "VARS" ) vars <- TRUE

    if (s == "TDEF") nt = as.integer(ctl[[i]][which(ctl[[i]]!="")][2])

  }
  
  # Use time steps from file if n is not given
  if (is.null(n) & !is.null(nt)) n = nt 
  time = c(1:n)

  # Determine total number of variables
  nvar <- length(nms)
  
  ## Done loading variable names
  
  newdata <- file(file.path(fldr,file), "rb")
  datavals <- readBin(newdata, "double", size = 4, n = nvar*n, endian = "little")
  close(newdata)
  
  dat <- as.data.frame(t(array(datavals,dim=c(nvar,n))))
  names(dat) <- nms
  dat$time = time 

  cat("\n Loaded: ",file.path(fldr,file),"\n")
  
  return(dat)
}

## Load 1d/2d/3d binary file
#' @export
load_grads <- function(fldr="/scratch/01/andrey/PPSCM",ctl="grads/history.ctl",file="OUT/history.dat",
                        n=NULL,nx=NA,ny=NA,mynms=c("ma","zs","zb","b"))
{
  ## Load variable names from CTL file
  ctl <- scan(file.path(fldr,ctl),what="character",sep="\n")
  ctl <- strsplit(ctl,split=" ")
  nms <- vector(mode="character")
  vars <- FALSE
  nt = NULL
  for ( i in 1:length(ctl)) {
    s <- ctl[[i]][1]
    if ( s == "ENDVARS" )vars <- FALSE
    if ( vars ) nms <- c(nms,s)
    if ( s == "VARS" ) vars <- TRUE

    if (s == "TDEF") nt = as.integer(ctl[[i]][which(ctl[[i]]!="")][2])
  }
  
  # Use time steps from file if n is not given
  if (is.null(n) & !is.null(nt)) n = nt 
  time = c(1:n)

  # Determine total number of variables and data fields
  nvar <- length(nms)
  nt   <- n
  ntot <- prod(c(nx,ny,nvar),na.rm=T)

  ## Done loading variable names
  
  newdata <- file(file.path(fldr,file), "rb")
  datavals <- readBin(newdata, "double", size = 4, n = nt*ntot, endian = "little")
  close(newdata)
  
  if ( is.na(nx) ) {  # 0D (time series)
  
    dat <- as.data.frame(t(array(datavals,dim=c(nvar,nt))))
    names(dat) <- nms

  } else {    # Multi-dimensional!!
    
    dat   <- list()
    empty <- array(NA,dim=c(nt,nx,ny))
    
    tmp <- array(datavals,dim=c(nx,ny,nvar,nt))
    
    for ( q in 1:nvar) {
      nm <- nms[q]

      if ( nm %in% mynms ) {
        dat[[nm]] <- empty
      for (k in 1:nt) dat[[nm]][k,,] <- tmp[,,q,k]
      }
      
    }
    
    #names(dat) <- nms

    # Add lat/lon
    dat$lat <- seq(from=  21,by=0.75,length.out=ny)
    dat$lon <- seq(from=-180,by=1.5, length.out=nx)
    
  }

  dat$time = time 
    
  cat("\n Loaded: ",file.path(fldr,file),"\n")
  
  return(dat)
}

## Load two binary files (from climber-sicopolis)
get.binaries <- function(fldr="/scratch/01/andrey/PPSCM",n,t0=0,accel=2,x1=NA,twoD=FALSE,x2=NA)
{
  
  ndat  <- n*1000/accel
  ndat2 <- n*1000
  n2d   <- n

  time1 <- c(1:ndat)*1e-3*accel+t0
  time2 <- c(1:ndat2)*1e-3 + t0
  time2d <- c(1:n2d) + t0

  if (is.na(x1[1])) x1 <- time1
  
  datclim <- get.binary2(fldr=fldr,ctl="grads/history.ctl",file="OUT/history.dat",n=ndat)
  datice  <- get.binary2(fldr=fldr,ctl="grads/time.ctl",file="sico_out/0d.ser",n=ndat2)
  
  dat2d <- NA
  if (twoD) dat2d   <- get.binary3(fldr=fldr,ctl="grads/2d_ser.ctl",file="sico_out/2d.ser",n=n2d,nx=241,ny=87)

  # Filter
  datclim <- sim.filter(datclim,x1=x1,x0=time1)
  datice  <- sim.filter(datice,x1=x1,x0=time2)
  
  tmp <- dat2d
  x0 <- c(1:dim(dat2d$zs)[1]) + t0
  kk <- which(x0 %in% x2)

  nms <- names(dat2d)
  for ( q in 1:length(nms) ) {
    nm <- nms[q]
    if ( ! nm %in% c("lat","lon") ) {
      tmp[[nm]] <- dat2d[[nm]][kk,,]
    }
  }
  dat2d      <- tmp
  dat2d$time <- x0[kk]

  # Combine series
  datice$time <- NULL
  dat <- cbind(datclim,datice)
  
  return(list(series=dat,dat2d=dat2d))
}
