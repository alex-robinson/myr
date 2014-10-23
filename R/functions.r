
plot.blank <- function(mar=NA)   
{ ## plot an empty space, in which text can be added
  def.par <- par(no.readonly=TRUE)
  if (is.na(mar[1])) mar <- c(1,1,1,1)
  par(mar=mar)
  plot(c(0,1),c(0,1),type="n",axes=FALSE, xlab="", ylab="")
  
  par(def.par$mar) # return to normal
}

darker <- function(col,percent=10)
{
  c <- rep("#000000",length(col))
  for (i in 1:length(col)) {
    cc <- col2rgb(col[i])/255 + (percent/100)
    cc[cc>1] <- 1; cc[cc<0] <- 0
    c[i] <- rgb(t(cc))
  }
  
  return(c)
}

alpha <- function(col,percent=50)
{
  c <- rep("#000000",length(col))
  for (i in 1:length(col)) {
    cc <- col2rgb(col[i])/255
    c[i] <- rgb(t(cc),alpha=percent/100)
  }
  
  return(c)
}

colmix <- function(col,col0="grey80",percent=50) {
    newcol = col
    for (q in 1:length(col)) newcol[q] = colorRampPalette(c(col0,col[q]))(100)[percent]
    return(newcol)
}

shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
                       theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) 
{

  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')

  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  text(xy$x, xy$y, labels, col=col, ... )
}

myfigure <- function(fldr=".",file="Rplot",date=TRUE,type="pdf",engine="cairo",
                     width=NULL,height=NULL,units="mm",asp=1,pointsize=12,res=300,
                     cex=1,cex.lab=1,cex.axis=1,bg="white",onefile=TRUE)
{
  
    # Make filename
    file = paste(file,".",type,sep="")
    if (date == TRUE) file = paste(today,"_",file,sep="")
    file = file.path(fldr,file)

    host = system("hostname",intern=TRUE)
    os   = system("uname",intern=TRUE)

    # If running on a mac, make sure engine is quartz!
    if ( os == "Darwin" ) engine = "quartz"

    # Determine width/heights in inches
    if ( is.null(width) & is.null(height) ) {  # Use default height, determine win via asp
        width  = 189  # Default width for pointsize 12
        height = width/asp
    } else if ( is.null(height) ) {  # only width specified, determine height
        height = width/asp
    } else if ( is.null(width) ) {  # only height specified, determine width
        width = asp*height
    } else {                    # height and width specified, determine asp
        asp = width/height
    }

    # Convert quantities if input was not inches
    cat(type,":",file,"\n")
    cat("width=",width,", height=",height," (",units,".) \n",sep="")
    conv = 1.0
    if ( units == "mm" ) conv = 0.0393700787
    if ( units == "cm" ) conv = 0.393700787
    hin = height*conv
    win = width*conv 
    #cat("width =",win,", height =",hin," (in.)\n")

    if (FALSE & os %in% c("Darwin") & type %in% c("png","jpg","tiff","pdf","ps")) {
        
        cat("Quartz plotting","\n")
        quartz(file=file,type=type,width=win,height=hin,pointsize=pointsize,dpi=res)

    } else if ( type == "png" ) {

        cat("engine = ",engine,"\n")
        png(file,width=win,height=hin,units="in",pointsize=pointsize,res=res,type=engine)

    } else if ( type == "jpg" ) {

        jpeg(file,width=win,height=hin,units="in",pointsize=pointsize,res=res,type=engine)

    } else if ( type == "tiff" ) {

        tiff(file,width=win,height=hin,units="in",pointsize=pointsize,res=res,type=engine)

    } else if ( type == "pdf" ) {

        if (engine %in% c("cairo","cairo1")) {
            cat("**cairo_pdf","\n")
            cairo_pdf(file,width=win,height=hin,pointsize=pointsize,onefile=onefile)
        } else {
            pdf(file,width=win,height=hin,pointsize=pointsize,onefile=onefile)
        }

    } else if ( type == "svg" ) {

        svg(file,width=win,height=hin,pointsize=pointsize)

    } else if ( type == "fig" ) {

        xfig(file,width=win,height=hin,pointsize=pointsize)

    } else {

        cairo_ps(file,width=win,height=hin,pointsize=pointsize)

    }

    par(bg=bg,cex=cex,cex.axis=cex.axis,cex.lab=cex.lab,tcl=0.2,mgp=c(2.5,0.3,0),las=1)

    return(win)
}


mylegend <- function(breaks,col,units="mm",x=c(0,1),y=c(0,1),at=NULL,labels=NULL,
                     xlab="",ylab="",xlim=NULL,ylim=NULL,zlim=range(breaks),
                     cex=1,cex.lab=1,new=TRUE,vertical=TRUE,line=1.8,
                     asp=1,mgp=c(3,0.5,0),col.axis="grey10",...)
{
    n      = length(breaks)    
    ynorm  = (breaks - min(breaks))
    ynorm  = ynorm / max(ynorm)
    y00    = ynorm[1:(n-1)]
    y11    = ynorm[2:n]
    x00    = rep(0,n)
    x11    = rep(1,n) 

    if ( vertical ) {
      x0   = x00
      x1   = x11 
      y0   = y00
      y1   = y11 
      xlim = c(0,1)
      ylim = zlim 
      ax   = 4
    } else {
      x0   = y00
      x1   = y11
      y0   = x00
      y1   = x11
      xlim = zlim
      ylim = c(0,1)
      ax   = 1
    }

    xlim0 = range(x0,x1)
    ylim0 = range(y0,y1)

    par(new=new,xpd=NA,xaxs="i",yaxs="i",...)
    plot( xlim0,ylim0, type="n",axes=F,ann=F,cex=cex)
    rect(x0,y0,x1,y1,col=col,border=col,lwd=1)

    par(new=TRUE,xpd=NA,xaxs="i",yaxs="i",...)
    plot(xlim,ylim,type="n",axes=F,ann=F,cex=cex)
    axis(ax,at=at,labels=labels,mgp=mgp,tcl=-0.1,col=col.axis,col.axis=col.axis,cex.axis=cex)
    box(col="grey10")

    mtext(side=1,line=line,xlab,cex=cex.lab)
    mtext(side=2,line=line,ylab,cex=cex.lab)

    par(xpd=FALSE)
}

my.par  <- function(mar=c(3.2,3.3,1,1),xaxs="i",yaxs="i",tcl=0.4,mgp=c(2.5,0.3,0),las=1,...)
{
  par(...,mar=mar,tcl=tcl,mgp=mgp,las=las,xaxs=xaxs,yaxs=yaxs)
}

my.axis <- function(side=1,at=NULL,tcl=0.4,mgp=c(2.5,0.25,0),minticks=2,grid=FALSE,...)
{
  
  if ( side %in% c(2,4)) mgp[2]=mgp[2]*1.3
  par(mgp=mgp,las=1)
  #if (side %in% c(1,3)) par(mgp=mgph)
  #if (side %in% c(2,4)) par(mgp=mgpv)
  
  # Calculate the main axis
  at.main = at
  if (is.null(at.main)) {
    at.main = axis(1,labels=F,tick=F)
  }

  # If desired, add ticks between standard values
  if (minticks > 1) {
    dx1 = diff(at.main)[1]/minticks
    at.min = axis(side=side,at=seq(range(at.main)[1]-10*dx1,range(at.main)[2]+10*dx1,by=dx1),labels=FALSE,tick=FALSE)
  } else {
    at.min = at.main
  }
  
  if (grid & side %in% c(1,3)) abline(v=at.min,lwd=1,lty=3,col=8)
  if (grid & side %in% c(2,4)) abline(h=at.min,lwd=1,lty=3,col=8)
  
  # Now actually draw ticks and labels...
  if (!grid & minticks > 1)
    axis(side=side,at=at.min, tcl=tcl/2,labels=FALSE)
  
  # Now plot main axes with labels, etc...
  axis(side=side,at=at.main,tcl=tcl,...)

  return(at.min)
}

my.plot <- function(...,axes=c(1,2,3,4),box=TRUE,grid=TRUE,
                    xlab="",ylab="",xline=1.7,yline=1.9,
                    minticks=2)
{
  #my.par(las=1)
  plot(...,axes=F,ann=F)
  #if (grid) grid()

  if ( 1 %in% axes) my.axis(side=1,minticks=minticks,grid=grid)
  if ( 2 %in% axes) my.axis(side=2,minticks=minticks,grid=grid)
  if ( 3 %in% axes) my.axis(side=3,labels=F,minticks=minticks)
  if ( 4 %in% axes) my.axis(side=4,labels=F,minticks=minticks)

  if (box) box()
  
  if ( !xlab=="" ) title(xlab=xlab,line=xline)
  if ( !ylab=="" ) title(ylab=ylab,line=yline)

}

fliplat <- function(dat,nc)
{ # Reverse the latitude dimension of all relevant variables of a data set
  
  vnms = names(nc$var)

  nms = NULL
  for (q in 1:length(vnms)) { 
    vnm = vnms[q]
    nd = nc$var[[vnm]]$ndim 
    latnow = FALSE
    for ( d in 1:nd ) { 
      if ("lat" %in% nc$var[[vnm]]$dim[[d]]$name) 
            nms = c(nms,vnm)
    }
  }

  cat("fliplat:",nms,"\n")

  # Latitude should be a vector, reverse it as needed
  if ( dat$lat[1] > dat$lat[2] ) {

    jj = rev(c(1:length(dat$lat)))
    dat$lat = dat$lat[jj]
    
    # Now reverse all variables 
    for ( q in 1:length(nms) ) {
      nm = nms[q]
      dims = dim(dat[[nm]])
      
      if ( nm == "lat_bnds" ) {

        dat[[nm]] = dat[[nm]][,jj]

      } else if ( length(dims)== 1 ) {

        dat[[nm]] = dat[[nm]][jj]
      
      } else if ( length(dims)== 2 ) {
        
        # Longitude is the first dimension
        # Latitude is the second dimension
        dat[[nm]] = dat[[nm]][,jj]

      } else if ( length(dims)==3 ) {
        
        # Longitude is the first dimension
        # Latitude is the second dimension
        dat[[nm]] = dat[[nm]][,jj,]

      } else if ( length(dims)==4 ) {
        
        # Longitude is the first dimension
        # Latitude is the second dimension
        dat[[nm]] = dat[[nm]][,jj,,]

      }

    } 
  
  }
  
  return(dat)

}

lon360to180 <- function(lon)
{
  # Adjust longitude
  lonX = lon
  i = which(lon > 180)
  lonX[i] = lonX[i] - 360
  i1 = which(lonX<0)
  i0 = which(lonX>0)
  ii = c(i1,i0)
  
  return(lonX)
}

shiftlons <- function(lon,lon180=TRUE)
{
    if (lon180 & range(lon,na.rm=TRUE)[2] > 180) {

    # Store indices of points in xlim range
    ii0 = which(lon >= xlim[1] & lon <= xlim[2])

    # Longitude
    lonX = lon
    i = which(lon > 180)
    lonX[i] = lonX[i] - 360
    i1 = which(lonX< 0)
    i0 = which(lonX>=0)
    ii = c(i1,i0)
    lonX = lonX[ii]

  } else if (!lon180 & range(lon,na.rm=TRUE)[1]<0) {

    # Longitude
    i0 = which(lon < 0)
    i1 = which(lon >=0)
    ii = c(i1,i0) 
    lonX = lon[ii]
    i = which(lonX < 0)
    lonX[i] = lonX[i] + 360
    
    #cat("lonX ",lonX,"\n")
  }

  return(list(lon=lonX,ii=ii))
}

shiftlon <- function(dat,nc)
{ # Shift longitude from 0:360 => -180:180 for
  # all related variables 
  
  vnms = names(nc$var)

  nms = NULL
  for (q in 1:length(vnms)) { 
    vnm = vnms[q]
    nd = nc$var[[vnm]]$ndim 
    latnow = FALSE
    for ( d in 1:nd ) { 
      if ("lon" %in% nc$var[[vnm]]$dim[[d]]$name) 
            nms = c(nms,vnm)
    }
  }
  
  cat("shiftlon:",nms,"\n")

  if (max(dat$lon) > 180 ) {    # Shift from 0:360 to -180:180

    # Longitude should be a vector
    ii = c(1:length(dat$lon))

    # Longitude
    lonX = dat$lon
    i = which(dat$lon > 180)
    lonX[i] = lonX[i] - 360
    i1 = which(lonX<0)
    i0 = which(lonX>0)
    ii = c(i1,i0)
    
    dat$lon = lonX[ii]
  
    # Now reverse all related variables 
    for ( q in 1:length(nms) ) {
      nm = nms[q]
      dims = dim(dat[[nm]])

      if ( nm == "lon_bnds" ) {

        dat[[nm]] = dat[[nm]][,ii]

      } else if ( length(dims)== 1 ) {

        dat[[nm]] = dat[[nm]][ii]
      
      } else if ( length(dims)== 2 ) {
        
        # Longitude is the first dimension
        # Latitude is the second dimension
        dat[[nm]] = dat[[nm]][ii,]

      } else if ( length(dims)==3 ) {
        
        # Longitude is the first dimension
        # Latitude is the second dimension
        dat[[nm]] = dat[[nm]][ii,,]

      } else if ( length(dims)==4 ) {
        
        # Longitude is the first dimension
        # Latitude is the second dimension
        dat[[nm]] = dat[[nm]][ii,,,]

      }

    } 
    
  }

  return(dat)

}


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

rmse <- function(x,ii=c(1:length(x)),round=NA)
{
  # Filter for desired points
  x <- x[ii]
  
  ii1 <- which(!is.na(x))
  rmse <- sqrt(sum(x[ii1]^2)/length(x[ii1]))
  
  if (!is.na(round)) rmse <- round(rmse,round)
  
  return(rmse)
  
}

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

# For countries and maps
require(maptools)
data(wrld_simpl)
require(gpclib)
# require(rgeos)
gpclibPermit()
wrld360 = nowrapRecenter(wrld_simpl,avoidGEOS=TRUE)

world_countries <- function() return(as.vector(wrld_simpl$NAME))

mask_by_region <- function(lon,lat,region=NULL)
{   # http://gis.stackexchange.com/questions/75033/given-a-lat-and-lon-identify-if-the-point-is-over-land-or-ocean

    if (is.null(region)) region = expand.grid(lon=lon,lat=lat)

    #region = rbind(region,region[1,])

    # Expand the grid array dimensions (lon,lat) to a set of data points
    points <- expand.grid(lon=lon,lat=lat)

    # Determine which points in the grid are inside of the regional polygon
    ii <- which(point.in.polygon(points[,1],points[,2],region[,1],region[,2]) > 0)

    # Make a mask of points inside the region
    mask = array(0,dim=c(length(lon),length(lat)))
    mask[ii] = 1

    return(mask)
}

mask_by_country <- function(lon,lat,countries="United States")
{   # http://gis.stackexchange.com/questions/75033/given-a-lat-and-lon-identify-if-the-point-is-over-land-or-ocean

    lon[lon>=180] = lon[lon>=180]-360 

    ## Create a SpatialPoints object
    points <- expand.grid(lon, lat)  # Note that I reversed OP's ordering of lat/long
    pts <- SpatialPoints(points, proj4string=CRS(proj4string(wrld_simpl)))

    ## Find which points fall over land
    qq = wrld_simpl$NAME %in% countries 

    ii <- !is.na(over(pts, wrld_simpl[qq,])$FIPS)

    mask = array(0,dim=c(length(lon),length(lat)))
    mask[ii] = 1

    return(mask)
}

my.in.land.grid <- function(lon,lat)
{ 
  countries = world_countries()
  mask = mask_by_country(lon,lat,countries)

  return(mask)
}


## Regional indices on Earth
regional_indices <- function(lon,lat,mask=NULL,file_srex="SREX_regions_adj.txt")
{ # Input vectors of lon and lat, returns grids of pre-defined regions
  
  # Get land mask
  #if (is.null(mask)) mask = in.land.grid(list(x=lon,y=lat))
  mask = my.in.land.grid(lon,lat)

  # Define polygons of regions of interest
  regions = list()
  
  # Global and global limited to 70 deg
  regions$glob    = data.frame(x = c(-180.00, -180.00, 180.00, 180.00),
                               y = c( -90.00,   90.00,  90.00, -90.00) )
  regions$glob70  = data.frame(x = c(-180.00, -180.00, 180.00, 180.00),
                               y = c( -70.00,   70.00,  70.00, -70.00) )
  
  # Hemispheres (North, South)
  regions$NH      = data.frame(x = c(-180.00, -180.00, 180.00, 180.00),
                               y = c(   0.00,   90.00,  90.00,   0.00) )
  regions$SH      = data.frame(x = c(-180.00, -180.00, 180.00, 180.00),
                               y = c(   0.00,  -90.00, -90.00,   0.00) )
  regions$WH      = data.frame(x = c(-180.00,    0.00,   0.00,-180.00),
                               y = c( -90.00,  -90.00,  90.00,  90.00) )
  regions$EH      = data.frame(x = c(   0.00,  180.00, 180.00,   0.00),
                               y = c( -90.00,  -90.00,  90.00,  90.00) )
  
  # Latitudinal circles
  regions$PRN     = data.frame(x = c(-180.00, -180.00, 180.00, 180.00),
                               y = c(  65.00,   90.00,  90.00,  65.00) )
  regions$PRS     = data.frame(x = c(-180.00, -180.00, 180.00, 180.00),
                               y = c( -65.00,  -90.00, -90.00, -65.00) )

  # Americas (North and South), Europe, Africa, Asia, Australia and New Zealand
  regions$NAM     = data.frame(x = c(-168.97, -169.79, -84.07, -34.45, -63.16, -80.79),
                               y = c(  77.24,   49.95,   4.57,  36.94,  71.84,  77.56) )
  regions$SAM     = data.frame(x = c( -71.77,  -87.76, -85.30, -60.29, -26.66, -71.36),
                               y = c(  14.09,   -4.63, -58.89, -58.89,  -6.53,  15.36) )
  regions$EUR     = data.frame(x = c( -24.80,  3.13,  9.78, 26.32, 26.03, 31.05, 39.62, 51.15, 50.56, 65.66, 69.81, 69.81, 15.54, -28.94),
                              y = c(   35.95, 36.64, 39.49, 34.01, 39.37, 42.60, 42.60, 42.60, 53.17, 53.28, 76.88, 82.01, 82.24,  64.91) )
  regions$AS      = data.frame(x = c( 180.00, 67.98, 69.04, 65.85, 50.43, 50.43, 35.53, 25.43, 27.02, 33.94, 43.51, 58.40, 113.19, 135.00, 180.00),
                               y = c(  90.00, 90.00, 73.88, 53.50, 53.50, 41.42, 43.28, 39.93, 34.33, 27.99, 11.94, 12.69, -15.30,   4.10,  39.18) )
  regions$AFR     = data.frame(x = c(  -1.23, -10.25, -23.79, -15.99, 23.38, 56.60, 53.32, 43.88, 33.22, 26.25, 10.25),
                               y = c(  36.62,  35.99,  21.71,  -0.82,-38.26,-21.76, 12.82, 11.87, 29.01, 33.45, 39.80) )
  regions$ANZ     = data.frame(x = c( 130.21, 106.81, 108.40, 179.15, 178.62),
                               y = c(   1.12, -19.03, -53.73, -54.48,   4.10) )
  
  regions$MED     = data.frame(x = c( -11.00,  45.00, 45.00,-11.00),
                               y = c(  26.00,  26.00, 48.00, 48.00) )

  ## Load SREX regions (used in http://www.ipcc-wg2.gov/SREX/)
  tmp = read.table(file_srex,header=T)
  tmp$lon.s = lon360to180(tmp$lon.s)
  tmp$lon.e = lon360to180(tmp$lon.e)
  tmp$lon.e[ which(tmp$lon.s > 0 & tmp$lon.e < 0) ] = 180.0

  nn = formatC(c(1:dim(tmp)[1]), width = 2, format = "d", flag = "0")
  srexnms = paste("SREX",nn,sep="")
  for (q in 1:length(nn)) {
    nm = srexnms[q]
    regions[[nm]] = data.frame( x = c(tmp$lon.s[q],tmp$lon.s[q],tmp$lon.e[q],tmp$lon.e[q]),
                                y = c(tmp$lat.s[q],tmp$lat.e[q],tmp$lat.e[q],tmp$lat.s[q]) )
  }
  
  ## Refine these to 12 continental regions
  

  # Now get indices of each region (using in.poly.grid)
  inds = list()
  nms = names(regions)
  for (q in 1:length(regions)) {
    # inds[[q]] = in.poly.grid(list(x=lon,y=lat), regions[[q]])
    inds[[q]] = mask_by_region(lon,lat,regions[[q]])
    if ( nms[q] %in% c("NAM","SAM","EUR","AFR","AS","ANZ","MED") |
         length(grep("SREX",nms[q]))>0 ) 
                 inds[[q]] = inds[[q]] & mask == 1
  }
  names(inds) = nms
  
  # Additional special regions
  inds$land    = mask == 1
  inds$NHland  = inds$NH & mask == 1
  inds$SHland  = inds$SH & mask == 1
  inds$ocean   = mask == 0 
  inds$NPocean = inds$PRN  & inds$ocean
  inds$NTocean = !inds$PRN & inds$NH & inds$ocean
  inds$STocean = !inds$PRS & inds$SH & inds$ocean
  inds$SPocean = inds$PRS  & inds$ocean
  
  ## CHECK REGIONS
  # col = tim.colors(length(inds))
  # image(x=lon,y=lat,z=mask,col=c("grey90","white"))
  
  # # Big regions
  # for(q in c(1:6)   +0) image(x=lon,y=lat,z=inds[[q]],col=c(NA,alpha(col[q],50)),add=T)
  # # Continents
  # for(q in c(1:6)   +6) image(x=lon,y=lat,z=inds[[q]],col=c(NA,alpha(col[q],50)),add=T)
  # # SRES
  # for (q in c(1:26)+12) image(x=lon,y=lat,z=inds[[q]],col=c(NA,alpha(col[q],50)),add=T)
  
  ## World Bank 3 regions ##
  MNA = c('Algeria','Bahrain','Djibouti','Egypt',
        'Iran (Islamic Republic of)','Iraq','Israel','Jordan','Kuwait',
        'Lebanon','Libyan Arab Jamahiriya','Morocco','Oman','Qatar',
        'Saudi Arabia','Syrian Arab Republic',
        'Tunisia','Untied Arab Emirates','Palestine','Yemen')

  LAC = c('Antigua and Barbuda','Argentina','Belize','Bolivia','Brazil', 
        'Chile','Colombia','Costa Rica','Dominica','Dominican Republic', 
        'Ecuador','El Salvador','Grenada','Guatemala','Guyana','Haiti',
        'Honduras','Jamaica','Mexico','Nicaragua','Panama','Paraguay', 
        'Peru','Saint Kitts and Nevis','Saint Lucia','Saint Vincent and the Grenadines',
        'Suriname','Uruguay','Venezuela')

  WBAL = c('Albania','Bosnia and Herzegovina',#'Kosovo',
         'The former Yugoslav Republic of Macedonia','Montenegro','Serbia')
  CAS  = c('Kazakhstan','Kyrgyzstan','Tajikistan','Turkmenistan','Uzbekistan')
  RUS  = "Russia"
  ECA  = c(WBAL,CAS,RUS)

  MEDNA = c('Algeria','Egypt','Israel','Palestine','Tunisia','Lebanon',
          'Libyan Arab Jamahiriya','Morocco','Turkey','Cyprus','Greece',
          'Albania','Croatia','The former Yugoslav Republic of Macedonia',
          'Bosnia and Herzegovina','Montenegro','Serbia','Italy','Spain','Malta')

  inds$MNA   = mask_by_country(lon,lat,countries=MNA)
  inds$LAC   = mask_by_country(lon,lat,countries=LAC)
  inds$WBAL  = mask_by_country(lon,lat,countries=WBAL)
  inds$CAS   = mask_by_country(lon,lat,countries=CAS)
  inds$RUS   = mask_by_country(lon,lat,countries=RUS)
  inds$ECA   = mask_by_country(lon,lat,countries=ECA)
  inds$MEDNA = mask_by_country(lon,lat,countries=MEDNA)
  WB3nms = c("MNA","LAC","ECA","WBAL","CAS","RUS","MEDNA")

  # Reduce to regions of interest
   nms = c("glob","glob70","land","NHland","SHland","NAM","SAM","EUR","AS","AFR","ANZ",
           srexnms,WB3nms,"MED",
           "ocean","NPocean","NTocean","STocean","SPocean")
  #nms = c("glob","land","NHland","ocean")
  inds = inds[nms]

  return(inds)
}

## Earth grid weighting ##
## AREA of grid boxes
gridarea <- function(lon,lat,Re=6371,method="cos")
{
  
  nx <- length(lon)
  ny <- length(lat)
  
  # Convert to radians
  latr <- lat*pi/180
  lonr <- lon*pi/180

  # Simple weighting based on cos(lat)...
  if (method == "cos") {

    a <- cos(latr)
    area <- matrix(rep(a,nx),byrow=TRUE,ncol=ny,nrow=nx)

  } else {
    
    # Take from: http://map.nasa.gov/GEOS_CHEM_f90toHTML/html_code/src/grid_mod.f.html
    #            2*PI*Re^2    {                                     }
    #    Area = ----------- * { sin( YEDGE[J+1] ) - sin( YEDGE[J] ) }
    #              IIGLOB     {                                     }
    # Re = radius of earth
    # YEDGE = position of latitude grid box edge
    # IIGLOB = number of longitudes
    
    dx <- (lonr[2] - lonr[1])/2
    dy <- (latr[1] - latr[2])/2
    
    area <- array(NA,dim=c(nx,ny))
    
    for ( j in 1:ny ) {
      a = (2*pi*Re^2/nx) * ( sin( latr[j]+dy ) - sin( latr[j]-dy ) )
      area[,j] = rep(a,nx)
    }

  }
  
  # Normalize the area to sum to 1
  area = area / base::sum(area)

  return(area)
}

mean.areawt <- function(var,area,ii=c(1:length(var)),mask=array(TRUE,dim=dim(var)),na.rm=T,...)
{

  # Limit area to masked area
  area[!mask] = NA

  # Reduce data to subset
  var  = var[ii]
  area = area[ii]

  # Remove NAs
  if (na.rm) {
    ii   = !is.na(var + area)
    var  = var[ii]
    area = area[ii]
  }
  
  ave = NA
  if (length(var) > 0) ave = sum(area * var)/sum(area)

  return(ave)
}

sd.areawt <- function(var,area,ii=c(1:length(var)),na.rm=T,normwt=T,...)
{ # Adapted from Hmisc package
  
  # Reduce data to subset
  var  = var[ii]
  area = area[ii]
  
  # std = wtd.var(x=as.numeric(var),weights=as.numeric(area),normwt=normwt)
  # std = sqrt(std)
  
  # Remove NAs
  if (na.rm) {
    ii   = !is.na(var + area)
    var  = var[ii]
    area = area[ii]
  }

  std = NA
  if (length(var)>0) {
    if (normwt) area = area * length(var)/sum(area)
    xbar = sum(area * var)/sum(area)
    std  = sum(area * ((var - xbar)^2))/(sum(area) - 1)
    std  = sqrt(std)
  }

  return(std)
}

normalize <- function(x,sd1=NA,ave=NA)
{
  if (is.na(sd1)) sd1 <- sd(x,na.rm=TRUE)
  if (is.na(ave)) ave <- mean(x,na.rm=TRUE)
  
  x1 <- (x-ave)/sd1
  return(x1)
}

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

findpoint <- function(dat,p=c(lat=37.82,lon=-25.73),ebs=3)
{ # Find the nearest grid point corresponding to a lat/lon location

  ii.point <- which( dat$lat >= (p[1]-ebs) &
                     dat$lat <= (p[1]+ebs) &
                     dat$lon >= (p[2]-ebs) &
                     dat$lon <= (p[2]+ebs) )
  
  d1 <- 1e6
  for (i in ii.point) {
    d <- calcdistw(p1=data.frame(lat=p[1],lon=p[2]),
                   p2=data.frame(lat=dat$lat[i],lon=dat$lon[i]))
    if ( d < d1 ) { d1 <- d; i1 <- i}
  }
  
  return(i1)
}

## LOADING BINARY (GRADS) FILES
# ----------------------------------------------------------------------
# Function: get.ctl
# Author  : Alex Robinson
# Purpose : Get info from a ctl file
# ----------------------------------------------------------------------
get.ctl <- function(fnm="clima.ctl")
{

  ctl <- system(paste("/home/robinson/python/sico/embfunctions.py",fnm),intern=TRUE)
  ctl <- strsplit(ctl," ")[[1]]

  gx = ctl[1]
  nx = as.integer(ctl[2])
  x0 = as.double(ctl[3])
  dx = as.double(ctl[4])
  ny = as.integer(ctl[5])
  y0 = as.double(ctl[6])
  dy = as.double(ctl[7])

  nvar = as.integer(ctl[8])
  vnms = ctl[9:length(ctl)]

  x = seq(from=x0,to=x0+(nx-1)*dx,by=dx)
  y = seq(from=y0,to=y0+(ny-1)*dy,by=dy)

  return(list(gx=gx,vnms=vnms,nvar=nvar,nxny=c(nx,ny),x=x,y=y))
}

# ----------------------------------------------------------------------
# Function: get.binary
# Author  : Alex Robinson
# Purpose : Get data from a binary output file
# ----------------------------------------------------------------------
get.binary <- function(fldr=".",gx="cntro.clima.gx",ctl="na",nx=151,ny=281,nk=1)
{

  if (ctl != "na") {
    ctl  <- get.ctl(file.path(fldr,ctl))

    gx   <- ctl$gx
    nvar <- ctl$nvar
    nx   <- ctl$nxny[1]
    ny   <- ctl$nxny[2]

    vnms <- ctl$vnms

  } else {
    nvar = length(vnms)
  }

  ndat = nx*ny
  newdata <- file(file.path(fldr,gx), "rb")
  datavals <- readBin(newdata, real(), size = 4, n = nvar*ndat*nk, endian = "little")
  close(newdata)

  cat("\n Loaded: ",file.path(fldr,gx),"\n")


  # Currently nk is not handled, also output should be in an array
  #var <- array(datvals,c(nx,ny,nk,nvar)

  var <- list()

  if (nk == 1) {

    # Obtain list of 2d variables
    for (q in 1:nvar)
    {
      i0 <- (q-1)*ndat + 1
      i1 <- i0 + ndat - 1
      var[[q]] <- matrix(datavals[i0:i1],nrow=nx,ncol=ny)
    }

    names(var) <- vnms

  } else {

    # Obtain list through time of list of 2d variables
    for (k in 1:nk) {

      var[[k]] <- list()

      for (q in 1:nvar)
      {
        i0 <- (k-1)*ndat*nvar + (q-1)*ndat + 1
        i1 <- i0 + ndat - 1
        var[[k]][[q]] <- matrix(datavals[i0:i1],nrow=nx,ncol=ny)
      }

      names(var[[k]]) <- vnms
    }

  }

  return(var)
}

## Load 1d binary file
get.binary2 <- function(fldr="/scratch/01/andrey/PPSCM",ctl="grads/history.ctl",file="OUT/history.dat",n=3)
{
  ## Load variable names from CTL file
  ctl <- scan(file.path(fldr,ctl),what="character",sep="\n")
  ctl <- strsplit(ctl,split=" ")
  nms <- vector(mode="character")
  vars <- FALSE
  for ( i in 1:length(ctl)) {
    s <- ctl[[i]][1]
    if ( s == "ENDVARS" )vars <- FALSE
    if ( vars ) nms <- c(nms,s)
    if ( s == "VARS" ) vars <- TRUE
  }
  
  # Determine total number of variables
  nvar <- length(nms)
  
  ## Done loading variable names
  
  newdata <- file(file.path(fldr,file), "rb")
  datavals <- readBin(newdata, real(), size = 4, n = nvar*n, endian = "little")
  close(newdata)
  
  dat <- as.data.frame(t(array(datavals,dim=c(nvar,n))))
  names(dat) <- nms
  
  cat("\n Loaded: ",file.path(fldr,file),"\n")
  
  return(dat)
}

## Load 1d/2d/3d binary file
get.binary3 <- function(fldr="/scratch/01/andrey/PPSCM",ctl="grads/history.ctl",file="OUT/history.dat",
                        n=3,nx=NA,ny=NA,mynms=c("ma","zs","zb","b"))
{
  ## Load variable names from CTL file
  ctl <- scan(file.path(fldr,ctl),what="character",sep="\n")
  ctl <- strsplit(ctl,split=" ")
  nms <- vector(mode="character")
  vars <- FALSE
  for ( i in 1:length(ctl)) {
    s <- ctl[[i]][1]
    if ( s == "ENDVARS" )vars <- FALSE
    if ( vars ) nms <- c(nms,s)
    if ( s == "VARS" ) vars <- TRUE
  }
  
  # Determine total number of variables and data fields
  nvar <- length(nms)
  nt   <- n
  ntot <- prod(c(nx,ny,nvar),na.rm=T)

  ## Done loading variable names
  
  newdata <- file(file.path(fldr,file), "rb")
  datavals <- readBin(newdata, real(), size = 4, n = nt*ntot, endian = "little")
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

## LOAD FILES WITH DATA FORMAT FROM ESRL ##
load_esrl <- function(filename,L=15)
{   # Make sure the first and last lines of the file are commented out first!

    tmp = read.table(filename)
    tmp[tmp==-99.99 | tmp==-99.9] = NA 
    dat = list(time=tmp[,1])
    dat$index.mon = as.array(t(tmp[,2:13]))
    dat$index = array(NA,dim=c(5,length(dat$time)))
    n = length(dat$time)
    tmp = dat$index.mon 
    tmp[12,2:n] = dat$index.mon[12,1:(n-1)]
    tmp[12,1]   = NA 
    dat$index[1,] = apply(tmp[c(12,1,2),],2,mean,na.rm=TRUE)
    dat$index[2,] = apply(tmp[c(3,4,5),],2,mean,na.rm=TRUE)
    dat$index[3,] = apply(tmp[c(6,7,8),],2,mean,na.rm=TRUE)
    dat$index[4,] = apply(tmp[c(9,10,11),],2,mean,na.rm=TRUE)
    dat$index[5,] = apply(tmp[c(1:12),],2,mean,na.rm=TRUE)

    dat$index.sm = dat$index 
    for (s in 1:5) dat$index.sm[s,] = ssatrend(dat$index[s,],L=L)

    return(dat)    
}

## Adding histogram to a plot
my.hist <- function(hi,freq=TRUE,b=FALSE,col='grey95',border='grey70',lwd=1,lty=1,filled="standard")
{
  
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

normalize_density <- function(dens,mids) 
{   # Returns the density renormalized (eg after multiplying probabilities together)

    dx = diff(mids[1:2])
    dens = dens / sum(dens) / dx
    return(dens)
}

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

join_pdfs1 <- function(dens2D,mids)
{   # Given density[mids,other], make a marginal pdf: density[mids]

    # Sum all pdfs together
    dens = apply(dens2D,1,sum)

    # Renormalize
    dens = normalize_density(dens,mids)

    out = list(x=mids,y=dens)
    return(out)
}

combine <- function(x,y,normalize=FALSE) 
{ # Function to combine all combinations of two vectors
  
  tmp <- expand.grid(x,y)
  out <- tmp[[1]]*tmp[[2]]
  
  if (normalize) out <- out/sum(out)
  
  cat("Combined values:",length(out),"\n")
  
  return(out)
}

plot.histweights <- function(dist)
{
  col <- rep(2,length(dist$x))
  col[dist$i95] <- "violet"
  col[dist$i66] <- 1
  col[dist$i10] <- 4
  plot(dist$x,dist$y,col=col,xlab="Temp (°C)",ylab="Weighting",xlim=c(0,6))

}
pointfit2 <- function(x0,y0,col=1,lwd=1,pch=1,lty=1)
{
  fit <- lm(y0~x0+I(x0^2))
  x1 <- sort(x0); y1 <- predict(fit,newdata=data.frame(x0=x1))
  points(x0,y0,pch=pch,col=col)
  lines(x1,y1,col=alpha(col,60),lwd=lwd,lty=lty)

}

## PLOTTING

# Default jet colors
jet.colors = c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F",
                 "yellow", "#FF7F00", "red", "#7F0000")

get.col = function(x,col=c("blue","red"),n=20,mid=NA,extend=10,ii=c(1:length(x)))
{
  nx = length(x)
  xlim = range(x,na.rm=T)
  
  if (!is.na(mid)) {
    # Caculate the xlim, so that mid is the midpoint! (eg, mid=0)
  }
  
  dxlim = diff(xlim)
  xlim[1] = xlim[1] - (extend/100)*dxlim
  xlim[2] = xlim[2] + (extend/100)*dxlim

  breaks = pretty(xlim,n)
  db = diff(breaks)[1]
  nb = length(breaks)

  palette = colorRampPalette(col)(nb-1)
  cols = rep(NA,nx)
  for ( i in 1:nx ) {
    j = which.min( abs(x[i] - breaks+db/2) )
    if ( length(j)==1 ) cols[i] = palette[j]
  }
  
  jj = which(! c(1:length(x) %in% ii))
  cols[jj] = NA
  return(list(breaks=breaks,palette=palette,col=cols))
}

