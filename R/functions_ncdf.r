

get.nc <- function(file,convert=1,vnms=NA,month=-1,missing=1e10,fliplat=TRUE,shiftlon=TRUE)
{
  nc <- open.ncdf(file)
  if (is.na(vnms[1])) vnms <- names(nc$var)
  dims = nc$dim 

  # Definet the initial list for output
  # including all the nc information
  out <- list(nc=nc)
  
  # Get dimensional variables, if they exist
  for ( q in 1:length(dims) ) {
    nm = dims[[q]]$name 
    out[[nm]]   <- get.var.ncdf(nc,nm)
  }

  # Load the variables
  for ( j in 1:length(vnms) ) {
    vnm <- vnms[j]
    
    v3      <- nc$var[[vnm]]
    varsize <- v3$varsize
    ndims   <- v3$ndims

    # Initialize start and count to read one timestep of the variable.
    start <- rep(1,ndims)   # begin with start=(1,1,1,...,1)
    start[ndims] <- 1       # change to start=(1,1,1,...,i) to read timestep i
    count <- varsize        # begin w/count=(nx,ny,nz,...,nt), reads entire var
    count[ndims] <- -1      # change to count=(nx,ny,nz,...,1) to read 1 tstep
    
    if ( month != -1 & v3$size[1] == 13 ) {
      start[1] <- month
      count[1] <- 1
    }
    
    #cat("vnm:",vnm,"\n")
    out[[vnm]] <- get.var.ncdf( nc, v3, start=start, count=count )
    out[[vnm]][abs(out[[vnm]]) >= missing] <- NA

  }

  # Make sure latest Vtot variable was obtained
  if (!is.null(out$Vtot)) {
  
    n <- length(out$Vtot)
    if ( is.na(out$Vtot[n]) ) out$Vtot[n] <- out$Vtot[n-1]
  
  }
  
  # Check if we are flipping a latitude dimension
  # (up to a three dimensional array assuming lat in second position!)
  if ( fliplat & "lat" %in% names(dims) )
    out = fliplat(dat=out,nc=nc)

  if ( shiftlon & "lon" %in% names(dims) )
    out = shiftlon(dat=out,nc=nc)

  cat("get.nc :",file,"\n")
  
  # Close files
  close.ncdf(nc)
  
  return(out)
}
