
#' @export
get_sector_contribs <- function(dat,sectors,nms=c("pp","snow","melt","runoff","smb"),dx=20)
{ ## Figure out sector contributions to SMB
  
  # Get conversion factor (mm=>Gt)
  conv <- (dx*1e3)^2*1e-12
  
  ns <- max(sectors) 
  ii0 <- which(dat$mask==0)
  
  sect <- as.data.frame(array(NA,dim=c(ns,length(nms)+2)))
  names(sect) <- c("sect","area",nms)
  
  for ( i in 1:ns ) {
    q <- which(sectors==i & dat$mask==0)
    
    sect$sect[i]   <- i
    sect$area[i]   <- length(q)/length(ii0)
    sect$pp[i]     <- sum(dat$pp[q])*conv
    sect$snow[i]   <- sum(dat$snow[q])*conv
    sect$melt[i]   <- sum(dat$melt[q])*conv
    sect$runoff[i] <- sum(dat$runoff[q])*conv
    sect$smb[i]    <- sum(dat$smb[q])*conv
  }
  
  return(sect)
}

#' @export
get_sector <- function(dat,latmid=72)
{ ## Return an array specifying the sector of each index
  
  zs  <- dat$zs
  lat <- dat$lat; lon <- dat$lon
  x   <- dat$xx;  y <- dat$yy
  
  ridge <- get_ridge(zs,x,y,lat,lon,ymin=-2900,ymax=-1700)
  
  ny <- dim(zs)[2]
  
  sector <- zs*NA
  
  for (q in 1:ny) {
    sector[,q] <- ifelse(lon[,q] > ridge$lon[q],1,2)
  }
  
  sector[lat > latmid & sector==2] <- 3
  sector[lat > latmid & sector==1] <- 4
  
  return(sector)
}

#' @export
get_ridge <- function(zs,x,y,lat,lon,ymin=-2900,ymax=-1700,spar=0.6)
{ ## Figure out where the central ridge line is based on maximum elevations
  ny <- dim(zs)[2]

  # Make some output vectors
  mid <- data.frame(j=numeric(ny),i=numeric(ny))
  mid$x <- mid$y <- mid$x2 <- mid$i2 <- NA
  mid$lat <- mid$lon <- NA
  
  for ( q in 1:ny ) {
    mid$y[q] <- y[1,q]; mid$j[q] <- q
    
    i <- which.max(zs[,q])
    mid$x[q] <- x[i,q]; mid$i[q] <- i
    
    i <- round(mean(order(zs[,q],decreasing=TRUE)[1:25]))
    mid$x2[q] <- x[i,q]; mid$i2[q] <- i
    
    mid$lat[q] <- lat[i,q]
    mid$lon[q] <- lon[i,q]
  }
  
  #r2 <- smooth.spline(mid$y,mid$x,spar=0.6)
  #mid$xsm <- r2$y; mid$ysm <- r2$x
  
  # Limit the ridge to certain y-values
  j0 <- which( round(mid$y,1) == ymin); j1 <- which( round(mid$y,1) == ymax)
  
  mid$x3 <- mid$x2; mid$y3 <- mid$y
  jj <- which(mid$y < ymin)
  mid$x3[jj] <- mid$x[j0];  mid$y3[jj] <- mid$y[j0]
  #mid$lon[jj] <- mid$lon[j0]; mid$lat[jj] <- mid$lat[j0]
  
  jj <- which(mid$y > ymax)
  mid$x3[jj] <- mid$x2[j1]; mid$y3[jj] <- mid$y[j1]
  #mid$lon[jj] <- mid$lon[j1]; mid$lat[jj] <- mid$lat[j1]
  
  return(mid)
}

#' @export
get_profile.dist <- function(dat,ii=which(dat$mask==0))
{ ## Figure out profiles of Greenland topography

  ii.water <- which(dat$mask==2)
  dist <- dat$mask*NA
  
  for ( i in ii ) {
    dists <- calcdist(data.frame(x=as.vector(dat$x[ii.water]),y=as.vector(dat$y[ii.water])),
                      data.frame(x=dat$x[i],y=dat$y[i]))
    dist[i] <- min(dists,na.rm=TRUE)
  }
  
  return(dist)
}

#' @export
get.slice <- function(dat,index=NA,time=NA,nms=c("mask","zs"),ebs=1e-3)
{ ## Get a specific time slice from a 2d dataset
  
  # Determine index of correct time
  if ( !is.na(time) ) {
    t <- which( abs(dat$time-time) < ebs)[1]
  } else {
    t <- index
    if (is.na(t)) t <- which.max(dat$time)
  }
  
  cat("get.slice: time=",dat$time[t],"\n")
  
  # First load generic grid variables (if they exist)
  out <- list()
  if ( "lat" %in% names(dat) ) out <- list(lat=dat$lat,lon=dat$lon,xx=dat$xx,yy=dat$yy)

  # Then load desired variables
  for ( q in 1:length(nms) ) {
    nm <- nms[q]
    if ( is.na(t) | nm %in% c("mask_hydro") ) {
      out[[nm]] <- dat[[nm]]
    } else {
      out[[nm]] <- dat[[nm]][,,t]
    }
  }
  
#   if ( "zs" %in% nms ) {
#     
#     i <- which.max(out$zs)
#     out$x.max <- grid$xx[i]; out$y.max <- grid$yy[i]
#     out$zs.max <- out$zs[i]
#     
#   }
  
  return(out)
}

#' @export
get_declineTiming <- function(time,Vol,V0=Vol[1],percent=c(20,50,80,85,90,100))
{ # Percent is the percent melted
  
  Vpp <- 100 * (Vol / V0)

  outV <- outt <- numeric(length(percent))*NA
  
  for ( q in 1:length(percent) ) {
    
    ii <- which( Vpp <= (100-percent[q]) )
    
#     if ( length(ii) == 0 ) ii <- c(1:length(Vpp))
    
    if ( length(ii) > 0 ) {
      i <- ii[which.max(Vpp[ii])]
      outV[q] <- Vpp[i]
      outt[q] <- time[i]
    }
    
  }
  
  out <- as.data.frame(as.list(c(outV,outt)))
  names(out) <- c(paste("V",percent,sep=""),paste("t",percent,sep=""))
  
  return(out)

}

#' @export
get_declineTiming2 <- function(time,Vol,V0=NULL,percent=20)
{ # Percent is the percent melted

  # Make sure Vol is a 2d array, where dim 1 is the simulation number
  # (for a vector, we only have one simulation)
  if (is.null(dim(Vol))) dim(Vol) = c(1,length(Vol))
  
  # How many simulations
  nsim = dim(Vol)[1]

  # Make sure the reference volume is defined  
  if (is.null(V0)) V0 = Vol[,1]
  
  # Get the volume percentage relative to reference volume
  Vpp <- 100 * (Vol / V0)
  
  out = list()
  out$timing  = array(NA,dim=c(nsim,length(percent)))
  out$volume  = array(NA,dim=c(nsim,length(percent)))
  out$percent = t(array(percent,dim=rev(dim(out$timing))))

  for (q in 1:nsim) {
    for (j in 1:length(percent)) {
      
      ii = which( Vpp[q,] <= (100-percent[j]) )

      if ( length(ii) > 0 ) {
        i = ii[which.max(Vpp[q,ii])]
        out$timing[q,j] = time[i]
        out$volume[q,j] = Vpp[q,i]
      }
    }
  }

  return(out)

}

#' @export
get.param.values <- function(fldrs,file=c("rembo.params","sico.params"),nms=NA,comment="#",skip=4)
{
  
  for ( q in 1:length(fldrs) ) {
    
    fldr <- fldrs[q]
    parameters1 <- matrix(scan(file.path(fldr,file[1]),what="character",comment.char=comment,skip=skip),byrow=TRUE,ncol=3)  
    
    parameters <- parameters1
    
    if ( length(file) > 1 ) {
      if ( file.exists(file.path(fldr,file[2])) ) {
        parameters2 <- matrix(scan(file.path(fldr,file[2]),what="character",comment.char=comment,skip=skip),byrow=TRUE,ncol=3)
        parameters  <- rbind(parameters1,parameters2)
      }
    }

    ## new ##
#     params <- as.list(parameters[,3])
#     names(params) <- parameters[,1]
#     for (i in 1:length(params)) {
#       p <- as.character(as.numeric(params[[i]]))
#       cat(names(params)[i]," : ","p =",p,"  : ","params[[i]] =",params[[i]],"\n")
#       if (!is.na(p) & p == params[[i]]) params[[i]] <- as.numeric(params[[i]])
#     }
    #######
    
#     params <- as.list((as.numeric(parameters[,3])))
    params <- as.list(parameters[,3])
    names(params) <- parameters[,1]
    
    nmstrings <- c("rcpname","ppsname","temper_file")
    
    ii <- which(! parameters[,1] %in% nmstrings )
    for (i in ii) {
      params[[i]] <- as.numeric(params[[i]])
    }
    
    if ( length( grep("PPS",params$temper_file) ) > 0 ) {
      params$ppsname <- substring(params$temper_file,first=14,last=18)
    }
    
    
    if (q == 1) { 
      all <- as.data.frame(params,stringsAsFactors=FALSE)
    } else {
      all <- rbind(all,as.data.frame(params,stringsAsFactors=FALSE))
    }
  
  }
  
  # Get indices of matching names
  pii <- c(1:dim(all)[2])
  if (!is.na(nms[1])) pii <- match(nms,names(all))
  
  bad <- which(is.na(pii))
  if (length(bad) > 0) {
    cat("Invalid parameter name(s):",nms[bad],"\n")
    pii <- pii[which(!is.na(pii))]
  }
  
  # Filter down to just the names I want
  all <- all[,pii]
  
  # Hack to account for changes to parameter names in rembo and sico
  if ( "margin_value" %in% names(all) ) all$margin_value <- NULL
  if ( "anf_frac" %in% names(all) )     all$anf_frac <- NULL
  if ( "rembo_start_file"  %in% names(all) ) { all$restart_file <- all$rembo_start_file; all$rembo_start_file <- NULL }
  
  # For rembo parameters only
  if ( "itm_c" %in% names(all) ) {
    if ( (!"itm_b" %in% names(all)) ) all$itm_b <- 0.0
    if ( (!"year_offset" %in% names(all)) ) all$year_offset <- 0.0
    if ( (!"T_diff" %in% names(all)) ) all$T_diff <- 0.0
  }
  
  # Modify some parameter names (backwards compatibility), if both rembo and sico used
  if ( "itm_c" %in% names(all) & "C_SLIDE_0" %in% names(all) ) {
    if ( "slide0" %in% names(all) ) { all$C_SLIDE_0 <- all$slide0; all$slide0 <- NULL }
    if ( "slide1" %in% names(all) ) { all$C_SLIDE_SEDI <- all$slide1; all$slide1 <- NULL }
    if ( "sea_level" %in% names(all) ) { all$SEA_LEVEL <- all$sea_level; all$sea_level <- NULL }
    if ( "ice_stream" %in% names(all) ) { all$ICE_STREAM <- all$ice_stream; all$ice_stream <- NULL }
    if ( "q_geo"  %in% names(all) ) { all$Q_GEO_0 <- all$q_geo; all$q_geo <- NULL }
    if ( "stream_cut.1" %in% names(all) ) { all$stream.cut.1 <- NULL }
  }
  
  #cat("fldr: ",fldrs,"\n",names(all),"\n")
  
  return(all)
  
}

#' @export
get.var <- function(file,ny=141,nx=76,what="double",comment="#")
{
  if ( what == "integer" ) {
    dat <- matrix(scan(file,comment=comment,what="character"),byrow=TRUE,nrow=ny,ncol=1)
    out <- matrix(NA,nrow=ny,ncol=nx)
    for ( q in 1:ny ) {
      out[q,] <- as.numeric(strsplit(dat[q,],"")[[1]])
    }
    out <- t(out)[,ny:1]
    
  } else {
    out <- t(matrix(scan(file,comment=comment),byrow=TRUE,nrow=ny,ncol=nx))[,ny:1]
  }
  
  return(out)
}

#' @export
topo.uplifted <- function(zs,zb,rho=910.0,rho_a=3300.0)
{
  return( zb + (rho/rho_a)*(zs-zb) )
}

#' @export
topo.stats <- function(zs,zb,mask=NA,dx=20,print=FALSE)
{
  
  # Define a mask
  if (is.na(mask[1])) {
    mask <- zs*0+2; mask[zs>=0] <- 1; mask[zs-zb>0] <- 0
  }
  
  conv.area <- 1e-6 * 1e-6  # m2 => 1e6 km2
  conv.vol  <- 1e-9 * 1e-6  # m3 => 1e6 km3
  
  # Get area in m2
  dx2 <- dx*dx*(1e3*1e3)
  
  H <- (zs - zb)
  ice.vol <- sum(H*dx2) * conv.vol
  
  tmp <- mask*0; tmp[mask==0] <- dx2
  ice.area <- sum(tmp) * conv.area
  
  tmp <- mask*0; tmp[mask<=1] <- dx2
  land.area <- sum(tmp) * conv.area
  
  if (print) {
    cat("\n","    ","  Land","  Ice","\n")
    cat("Area",round(land.area,3),round(ice.area,3),"\n")
    cat("Vol "," --- ",round(ice.vol,3),"\n")
    cat("\n","Elev. range:",round(min(zs),1),round(max(zs),1),"\n")
    cat("\n")
  }
  
  return(data.frame(ice.vol=ice.vol,ice.area=ice.area,land.area=land.area))
}

#' @export
stream.gen <- function(zs,fmax=1,zmax=500)
{ ## Generate a sliding mask, where the sliding factor
  ## increases with lower elevations
  
  # First make a mask of zeros (assuming no sliding anywhere)
  mask <- zs*0
  
  # Obtain points that are of interest, ie, less than zmax
  ii <- which(zs <= zmax)
  
  # Determine how strong the sliding factor is at each point
  # (the lower the elevation, the stronger the sliding factor
  mask[ii] <- fmax*(1-zs[ii]/zmax)
  
  # If the elevation is below zero, or the mask
  # erroneously gives a value greater than fmax, reset to fmax
  mask[zs <= 0 | mask > fmax] <- fmax
  
  return(mask)
}

#' @export
extend.hydrobasins <- function(grid,basins)
{  ## Calculate the hydrological basins for the whole grid
    
  # Store original data
  bas0  <- bas1 <- basins
  
  nx <- dim(grid$xx)[1]
  ny <- dim(grid$xx)[2]
  
  # Eliminate non-existent basins
  bas0[bas0==0] <- NA
  
  # Loop over all points and fill in basins
  # for points lacking a label (with the nearest neighbor)
  k <- 0
  n <- nx*ny
  
  for ( i in 1:nx ) {
    for (j in 1:ny ) {
      
      if (is.na(bas0[i,j])) {
        
        p0 <- list(x=grid$xx[i,j],y=grid$yy[i,j])
        p1 <- list(x=as.vector(grid$xx),y=as.vector(grid$yy))
        
        dists <- calcdist(p0,p1)
        
        qq <- order(dists)
        q <- qq[which(!is.na(bas0[qq]))][1]
        
        bas1[i,j] <- bas0[q]
      }
     
    }
    
    cat("Row", i,"of",nx,"\n")
  }

  # Return new basins
  return(bas1)

}

#' @export
monthly2daily.smooth <- function(Tm,day=c(1:360),method="linear")
{
  Tm <- as.numeric(Tm[tolower(month.abb)])
  d0 <- c(1:12)*30 - 15
  
  if ( method == "linear" ) {
    Td <- approx(d0,Tm,day,rule=2)
  } else if ( method == "spline" ) {
    Td <- spline(d0,Tm,xout=day,method="periodic")
  }
  
  return(Td$y)
}

#' @export
monthly2daily <- function(Tm,day=c(1:360))
{ ## Interpolate data in time: monthly => daily
  
  Tm <- as.numeric(Tm)
  
  ndm = 30
  nm  = 12
  
  n <- length(day)
  m0  = m1  = numeric(n)
  wt0 = wt1 = numeric(n)
  Td        = numeric(n)
  
  # Get length of arrays to fill, midpoint of month (in days)
  # and the total weight (total number of days in a month)
  mid = ndm / 2
   
  for (k in 1:n) {
    
    d <- day[k]
    
    for( m in 1:(nm+1)) {
      if ( m*ndm-mid > d ) break
    }
    m1[k] = m; m0[k] = m1[k]-1
    
    if ( m1[k] > 12 ) m1[k] =  1
    if ( m0[k] <  1 ) m0[k] = 12
    
#     wt1[k] = abs( mod(k-mid,ndm) )
    wt1[k] = abs( (d-mid) %% ndm )
    wt0[k] = ndm - wt1[k]
    
    wttot = wt0[k] + wt1[k]
    wt0[k] = wt0[k]/wttot; wt1[k] = wt1[k]/wttot
    
    Td[k] = wt0[k]*Tm[m0[k]] + wt1[k]*Tm[m1[k]]
  }
   
  return(Td)
  
}

#' @export
effectiveT <- function(T,sigma=5)
{
  inv_sqrt2   = 1.0/sqrt(2.0)
  inv_sqrt2pi = 1.0/sqrt(2.0*pi)

  inv_sigma   = 1.0/sigma

  Teff = sigma*inv_sqrt2pi*exp(-0.5*(T*inv_sigma)^2) +
               T*0.5*erfcc(-T*inv_sigma*inv_sqrt2)

  return(Teff)
}

#' @export
erfcc <- function(x)
{
  z = abs(x)
  t = 1.0/(1.0+0.5*z)

  erfc = t*exp(-z*z-1.265512230+t*(1.000023680+t*(0.374091960 +
    t*(0.096784180+t*(-0.186288060+t*(0.278868070+t*
    (-1.135203980+t*(1.488515870+t*(-0.822152230+
    t*0.170872770)))))))))

  #if (x < 0.0) { erfcc = 2.0-erfcc }
  erfc[x < 0.0] = 2.0-erfc[x < 0.0]

  return(erfc)
}

### RCM data loading

#' @export
dma2cma <- function(dma,area=1.74,conv=1e-2*area*1e6*1e-6)
{ # dma given in percent extent
  # convert percent to fraction
  # convert area to km2 from 1e6 km2
  # convert final cma from km2 to 1e6 km2
  return(sum(dma)*conv)
}

#' @export
get_cma_series <- function(dat,nm="MAR_MELT1B",mnm="MSK",mice=1,t0=1958,dx=25,nd=365,area=1.74)
{

  mask <- dat[[mnm]]; ii <- which(mask!=mice)
  n <- dim(dat[[nm]])[3]
  
  t <- seq(from=t0,by=1,length.out=n)
  
  ma <- array(0,dim=dim(dat[[nm]]))
  cma <- ama <- numeric(n)
  
  # Loop over each year to output the filtered mask and the CMA
  for ( j in 1:n ) {
    
    melt <- dat[[nm]][,,j]; melt[ii] <- 0
    
    # Store the array of melt area for this year
    ma[,,j] <- melt
    
    # Now determine CMA for this year
    cma[j] <- sum(melt)*(dx)^2 * 1e-6
    
    ama[j] <- 100* (cma[j] / nd) / area
  }
  
  return(list(time=t,cma=cma,ama=ama,ma=ma))
}

#' @export
Gt2sl <- 1/362   # 1mm = 362 Gt

#' @export
aline <- function(time=NA,var=NA,t0=1979,tf=2050,y0=0.7556875,slope=36e4*1e-6,sd=17e4*1e-6/3)
{
  # Get the fit here
  if (!is.na(time[1])) {
    ii <- which(time >= t0 & time <= tf)
    x00 <- time[ii]; y00 <- var[ii]
    fit <- lm(y00~x00)
    
    sd <- sd(fit$residuals)
    y0 <- fit$coefficients[1]; slope <- fit$coefficients[2]
    y0 <- y0 + t0*slope
  }
  # fit obtained
  
  # Get fit and confidence intervals
  t <- seq(from=t0,to=tf,by=1)
  y <- predict(fit,data.frame(x00=t),interval="confidence")
  y1 <- y[,2]; y2 <- y[,3]
  y <- y[,1]
  
#   y <- y0 + (t-t0)*slope
#   
#   y1 <- y0 + (t-t0)*(slope-sd)
#   y2 <- y0 + (t-t0)*(slope+sd)
  
  return(list(slope=slope,sd=sd,time=t,y=y,y1=y1,y2=y2))
}

#' @export
load.rcmdata <- function(file="../rcms/RCMGCM_mar.era40.txt",
                             time=NA,trange=c(1900,2500),trange.base=c(1958,2001))
{ ## Load temperature time-series and add to datasets
  
  # Read the file and extract the average of the months of interest
  tmp <- read.table(file,header=TRUE)
  names(tmp)[2:13] <- toupper(month.abb)
  tmp$tjja <- apply(tmp[,c("JUN","JUL","AUG")],FUN=mean,MARGIN=c(1))
  
  n <- length(tmp$DEC)
  tmp$DECm1 <- tmp$DEC
  tmp$DECm1[2:n] <- tmp$DEC[1:(n-1)]
  
  tmp$tdjf  <- apply(tmp[c("DECm1","JAN","FEB")],FUN=mean,MARGIN=1)
  tmp$tdjfm <- apply(tmp[c("DECm1","JAN","FEB","MAR")],FUN=mean,MARGIN=1)
  
  # Filter to the time range of interest
  ii <- which(tmp$time >= trange[1] & tmp$time <= trange[2])
  if (!is.na(time[1])) ii <- which(tmp$time %in% time)
  
  out <- as.list(tmp[ii,])
  
  # Determine the base temperature during the base time period and
  # subtract to get anomalies
  ii <- which(out$time >= trange.base[1] & out$time <= trange.base[2])
  
  ## WRONG - I needed to add ii here (18.05.2011)
  out$tjja0 <- mean(out$tjja[ii],na.rm=TRUE)
  out$dtjja <- out$tjja - out$tjja0
  
  out$tdjf0 <- mean(out$tdjf[ii],na.rm=TRUE)
  out$dtdjf <- out$tdjf - out$tdjf0
  
  if ( "snow" %in% names(out) ) out$precip <- out$snow + out$rain
  
  return(out)
}


## Load and process RCM data
#' @export
save.rcmdata <- function()
{
  
  # Load RCM daily data
  dmas <- read.table("../rcms/Fettweis_dma_1958-2001.txt",header=TRUE)
  dmas$Day <- dmas$Day*360/365   # Scale days to 360-day year
  
  dmas2 <- read.table("../rcms/Fettweis_dma_1979-2009.txt",header=TRUE)
  dmas2$Day <- dmas2$Day*360/365   # Scale days to 360-day year_offset
  
  ## Get CMA from daily dataset
  cma.thresh <- c(7.75,8.50,9.25)
  cma.mar <- c(dma2cma(dmas$MAR1a),dma2cma(dmas$MAR1b),dma2cma(dmas$MAR1c))
  cma.rac <- c(dma2cma(dmas$RACMO1a),dma2cma(dmas$RACMO1b),dma2cma(dmas$RACMO1c))
  
  nd <- 360; conv <- 100/1.74
  ama.mar <- conv* cma.mar / nd
  ama.rac <- conv* cma.rac / nd
  
  # Load yearly melt data (2D)
  dat <- get.nc("../rcms/MAR-RACMO2-SAT_melt_extent.nc")
  
  # MAR
  mar1b <- load.rcmdata(file="../rcms/RCMGCM_mar.era40.txt",trange.base=c(1958,2001))                         
  tmp   <- get_cma_series(dat,nm="MAR_MELT1B",mnm="MSK",mice=1,t0=1958,dx=25) 
  ii    <- which(tmp$time %in% mar1b$time)
  ii0 <- which(mar1b$time %in% tmp$time[ii])
  mar1b$cma <- mar1b$ama <- numeric(length(mar1b$time))*NA
  mar1b$ma <- array(NA,dim=c(dim(tmp$ma)[1:2],length(mar1b$time)))
  mar1b$cma[ii0]  <- tmp$cma[ii]
  mar1b$ama[ii0]  <- tmp$ama[ii]
  mar1b$ma[,,ii0] <- tmp$ma[,,ii]
  
  # RACMO
  rac1b <- load.rcmdata(file="../rcms/RCMGCM_racmo.era40.txt",trange.base=c(1958,2001))  
  tmp   <- get_cma_series(dat,nm="RAC_MELT1B",mnm="MSK",mice=1,t0=1958,dx=25)
  ii    <- which(tmp$time %in% rac1b$time)
  ii0 <- which(rac1b$time %in% tmp$time[ii])
  rac1b$cma <- rac1b$ama <- numeric(length(rac1b$time))*NA
  rac1b$ma <- array(NA,dim=c(dim(tmp$ma)[1:2],length(rac1b$time)))
  rac1b$cma[ii0]  <- tmp$cma[ii]
  rac1b$ama[ii0]  <- tmp$ama[ii]
  rac1b$ma[,,ii0] <- tmp$ma[,,ii]
  
  # SAT
  sat1b <- get_cma_series(dat,nm="SAT_MELT1B",mnm="MSK",mice=1,t0=1958,dx=25)
  # Make an artificial temperature series for satellite data
  ii1 <- which(mar1b$time %in% sat1b$time)
  ii2 <- which(rac1b$time %in% sat1b$time)
  
  nms <- c(toupper(month.abb),"tjja","tdjf","tdjfm")
  for (q in 1:length(nms)) {
    nm <- nms[q]
    sat1b[[nm]] <- apply( rbind(mar1b[[nm]][ii1],rac1b[[nm]][ii2]),FUN=mean,MARGIN=c(2) )
  }
  sat1b$tjja0 <- mean( sat1b$tjja[ which(sat1b$time >= 1958 & sat1b$time <= 2001) ] )
  sat1b$dtjja <- sat1b$tjja - sat1b$tjja0
  
  # Get base tjja for 1980-1999 like GCMs
  t0 <- 1980; t1 <- 1999
  ii <- which(mar1b$time >= t0 & mar1b$time <= t1)
  mar1b$tjja0gcm <- mean(mar1b$tjja[ii],na.rm=TRUE)
  ii <- which(rac1b$time >= t0 & rac1b$time <= t1)
  rac1b$tjja0gcm <- mean(rac1b$tjja[ii],na.rm=TRUE)
  ii <- which(sat1b$time >= t0 & sat1b$time <= t1)
  sat1b$tjja0gcm <- mean(sat1b$tjja[ii])
  dtjja0 <- mean(sat1b$dtjja[ii],na.rm=TRUE)
  
  # Add dummy field to mimic gcm fields
  mar1b$dtjja0 <- mar1b$dtjja
  rac1b$dtjja0 <- rac1b$dtjja
  sat1b$dtjja0 <- sat1b$dtjja
  
  ## Load GCM cmas and temperatures...
  
  # CMA data
  datcma <- read.table("../rcms/RCMGCM_proj.cma.txt",header=TRUE)
  datcma[,2:dim(datcma)[2]] <- datcma[,2:dim(datcma)[2]]*1e1
  nd.gcm <- 365
  
  # MAR - ECHAM5 - A1B
  MAR.echam5.A1B <- load.rcmdata(file="../rcms/RCMGCM_mar.echam5.A1B.txt",trange.base=c(1980,1999)) 
  ii <- which(datcma$time %in% MAR.echam5.A1B$time)
  MAR.echam5.A1B$cma <- datcma$MAR.echam5.A1B[ii]                      
  MAR.echam5.A1B$ama <- conv* MAR.echam5.A1B$cma / nd.gcm
  
  # Offset the anomaly dtjja to account for bias in GCM
  ii <- which(MAR.echam5.A1B$time >= t0 & MAR.echam5.A1B$time <= t1)
  dtjja1 <- mean(MAR.echam5.A1B$dtjja[ii])
  MAR.echam5.A1B$dtjja.offset <- dtjja1-dtjja0
  MAR.echam5.A1B$dtjja0 <- MAR.echam5.A1B$dtjja
  MAR.echam5.A1B$dtjja <- MAR.echam5.A1B$dtjja0 - MAR.echam5.A1B$dtjja.offset
  
  
  # MAR - ECHAM5 - E1
  MAR.echam5.E1 <- load.rcmdata(file="../rcms/RCMGCM_mar.echam5.E1.txt",trange.base=c(1980,1999)) 
  ii <- which(datcma$time %in% MAR.echam5.E1$time)
  MAR.echam5.E1$cma <- datcma$MAR.echam5.E1[ii]
  MAR.echam5.E1$ama <- conv* MAR.echam5.E1$cma / nd.gcm
  
  # Offset the anomaly dtjja to account for bias in GCM
  ii <- which(MAR.echam5.E1$time >= t0 & MAR.echam5.E1$time <= t1)
  dtjja1 <- mean(MAR.echam5.E1$dtjja[ii])
  MAR.echam5.E1$dtjja.offset <- dtjja1-dtjja0
  MAR.echam5.E1$dtjja0 <- MAR.echam5.E1$dtjja
  MAR.echam5.E1$dtjja <- MAR.echam5.E1$dtjja0 - MAR.echam5.E1$dtjja.offset
  
  
  # MAR - HADCM3 - A1B   
  MAR.hadcm3.A1B <- load.rcmdata(file="../rcms/RCMGCM_mar.hadcm3.A1B.txt",trange.base=c(1980,1999)) 
  ii <- which(datcma$time %in% MAR.hadcm3.A1B$time)
  MAR.hadcm3.A1B$cma <- datcma$MAR.hadcm3.A1B[ii]
  MAR.hadcm3.A1B$ama <- conv* MAR.hadcm3.A1B$cma / nd.gcm
  
  # Offset the anomaly dtjja to account for bias in GCM
  ii <- which(MAR.hadcm3.A1B$time >= t0 & MAR.hadcm3.A1B$time <= t1)
  dtjja1 <- mean(MAR.hadcm3.A1B$dtjja[ii])
  MAR.hadcm3.A1B$dtjja.offset <- dtjja1-dtjja0
  MAR.hadcm3.A1B$dtjja0 <- MAR.hadcm3.A1B$dtjja
  MAR.hadcm3.A1B$dtjja <- MAR.hadcm3.A1B$dtjja0 - MAR.hadcm3.A1B$dtjja.offset
   
  # Calculate projections
  pmar1b <- aline(time=mar1b$time,var=mar1b$cma,y0=NA,t0=1979,tf=2050)
  psat1b <- aline(time=sat1b$time,var=sat1b$cma,y0=NA,t0=1979,tf=2050)
  prac1b <- aline(time=rac1b$time,var=rac1b$cma,y0=NA,t0=1979,tf=2050)
  
  save(dmas, dmas2, cma.thresh, cma.mar, cma.rac, ama.mar, ama.rac,
       mar1b, rac1b, sat1b, MAR.echam5.A1B, MAR.echam5.E1,MAR.hadcm3.A1B,
       file="rcmdata.Rdata")
  
}

## NEW ##
#' @export
load_MAR0 <- function(fldr="data/rcms/MARv2/MARv2_ERA-40_1958-1999")
{
  
  vnms = c("ice","me","melt","rf","ru","sf","stt","su","tt")
  files = list.files(fldr)
  
  if (length(vnms) != length(files) ) {
    cat("Wrong files. Fix function 'load_MAR'!\n")
    cat(vnms,"\n")
    cat(files,"\n")
  }
  
  dat = list()
  for (q in 1:length(vnms)) {
    vnm = vnms[q]
    q0 = grep(paste("-",vnm,".dat",sep=""),files)

    if ( length(q0) > 0 ) {

      tmp = read.table(file.path(fldr,files[q0]))

      if (q==1) dat$time = tmp[[1]]
      dat[[vnm]] = tmp[,c(2:13)]
      names(dat[[vnm]]) = tolower(month.abb)
    }

  }
  
  cat("Loaded MAR folder:",fldr,"\n")

  return(dat)
}

#' @export
load_MAR <- function(fldr,time=c(1950:2100),time.base=c(1980:2000))
{
  icearea = 1744980  # km2

  # Now load MAR data set(s)
  fldr1 = NULL
  if (length(grep("rcp",fldr))>0) {
    fldr1 = fldr
    fldr1 = gsub("rcp26_2006-2100","histo_1965-2005",fldr1)
    fldr1 = gsub("rcp45_2006-2100","histo_1965-2005",fldr1)
    fldr1 = gsub("rcp85_2006-2100","histo_1965-2005",fldr1)
  }
  if (length(grep("ERA-40",fldr))>0) {
    fldr1 = fldr
    fldr1 = gsub("40_1958-1999","INTERIM_1979-2011",fldr1)
  }
  
  dat0 = list(load_MAR0(fldr))
  if (!is.null(fldr1)) dat0[[2]] = load_MAR0(fldr1)
  
  # dat = list(time=time,month=c(1:12),tjja=NA,dtjja=NA,cma=NA)
  dat = list(time=time,month=c(1:12))
  nt = length(time)
  dat$tt   = as.data.frame(array(NA,dim=c(nt,12)))
  dat$melt = as.data.frame(array(NA,dim=c(nt,12)))
  dat$smb  = as.data.frame(array(NA,dim=c(nt,12)))

  for (q in 1:length(dat0)) {
    
    tmp = dat0[[q]]
    
    # Filter to the time range of interest
    ii0 <- which(tmp$time %in% dat$time)
    ii1 <- which(dat$time %in% tmp$time)

    dat$tt[ii1,]   = tmp$tt[ii0,]
    dat$melt[ii1,] = tmp$melt[ii0,]
    dat$smb[ii1,]  = tmp$sf[ii0,] + tmp$rf[ii0,] - tmp$ru[ii0,] - tmp$su[ii0,]

  }
  
  dat$tjja = apply(dat$tt[,c(6,7,8)],FUN=mean,MARGIN=c(1))
  dat$cma  = apply(dat$melt*icearea,FUN=sum,MARGIN=c(1))
  dat$cma[dat$cma==0] = NA

  ## Get temp anomalies relative to base period
  ii = which(dat$time %in% time.base)
  T0 = mean(dat$tjja[ii])
  dat$dtjja = dat$tjja - T0
  
  ## Get REMBO base period if ERA-40
  if (length(grep("ERA-40",fldr))>0) {
    ii = which(dat$time %in% c(1958:2001)) 
    T0 = mean(dat$tjja[ii])
    dat$dtjja.rembo = dat$tjja - T0
  }
  
  return(dat)
}

###

gen.styles <- function(nm,n=length(nm),col=rep(1,n),lty=rep(1,n),
                       lwd=rep(1,n),pch=rep(23,n),pch.lwd=rep(1,n),cex=rep(1,n),bg=rep(NA,n))
{
  styles <- data.frame(nm=nm,col=col,lty=lty,lwd=lwd,pch=pch,pch.lwd=pch.lwd,bg=bg,cex=cex)
  styles$lty     <- as.numeric(styles$lty)
  styles$lwd     <- as.numeric(styles$lwd)
  styles$pch     <- as.numeric(styles$pch)
  styles$pch.lwd <- as.numeric(styles$pch.lwd)
  styles$bg      <- as.numeric(styles$bg)
  styles$cex     <- as.numeric(styles$cex)
  return(styles)
}

land.colors  <- colorRampPalette(c("tan","brown"),bias=5)
water.colors <- colorRampPalette(c("skyblue3","lightblue"))
grey.colors  <- colorRampPalette(c("grey50","grey80"))
abc.labels <- c("a","b","c","d","e","f","g")

panel.contour.grl <- function(x,y,subscripts,...,datasets,
                              extra=extra,extra2=extra2,ptext=NA,plabel=NA,pobs=NULL,latlon=TRUE)
{
  # Determine which data set is being plotted
  pnow <- panel.number()
  
  # Get info for current data set
  land <- extra$land; mask <- extra$mask; elev <- extra$elev; points <- extra$points
  var  <- extra$var
  vnm  <- extra$vnms[ min(pnow,length(extra$vnms)) ]
  
  cat("var: ",vnm,"\n")
  
  # Get current data set
  pdat <- datasets[[ min(pnow,length(datasets)) ]]
  out  <- pdat[[vnm]]  
  ma   <- pdat$mask; zs  <- pdat$zs
  xx   <- pdat$xx; yy  <- pdat$yy
  ii0 <- c(1:length(out))
  
  cat("n=",length(out),"\n")  
  
  # Plot actual variable for subscripts desired 
  #panel.levelplot(x,y,subscripts=subscripts,...)
  ii1 <- which(ma %in% var$mask)
  panel.contourplot(xx, yy, out,subscripts=ii1,
                    at=var$at, col.regions=var$cols,
                    col=var$contour.col,lty=var$contour.lty,lwd=var$contour.lwd,
                    region = TRUE, contour = var$contour.plot, labels = FALSE)
  
  # Plot land if desired
  if (land$plot) {
    ii1 <- which(ma %in% land$mask)
    panel.contourplot(xx, yy, zs,subscripts=ii1,
                      at=land$at, col.regions=land$cols,
                      region = TRUE, contour = FALSE, labels = FALSE)
  }
  
  # Plot the ocean
  if ( ocean$plot ) {
    ii1 <- which(ma %in% ocean$mask)
    panel.contourplot(xx, yy, zs,subscripts=ii1,
                      at=ocean$at, col.regions=ocean$cols,
                      region = TRUE, contour = FALSE, labels = FALSE)
  }
  
  # Plot land/ice outline from mask
  if (mask$plot) {
    # First plot the desired mask contours
    panel.contourplot(xx, yy, ma, subscripts=ii0,lty=1,col=1,lwd=2,
                      at=mask$at,
                      region = FALSE, contour = TRUE, labels = FALSE)
  }
  
  # Add elevation contours (thicker lines for round numbers..)
  if (elev$plot) {
    ii1 <- which(ma %in% elev$mask)
    panel.contourplot(xx, yy, zs,subscripts=ii1,lty=2,col="grey20",lwd=0.5,
                      at=elev$at,region = FALSE, contour = TRUE, labels = FALSE)
  }
  
  # Add additional contour field (eg, data for comparison
  if (!is.null(extra2)) {
    # Plot the desired mask contours
    panel.contourplot(xx, yy, extra2$mask, subscripts=ii0,lty=extra2$lty,col=extra2$col,lwd=extra2$lwd,
                      at=extra2$at,
                      region = FALSE, contour = TRUE, labels = FALSE)
  }
  
  # Add observational or other points of interest to plot
  if (points$plot) {
    #c2 <- rgb(t(1-col2rgb(points$col)))
    panel.points(points$x,points$y,col=points$col,pch=points$pch,fill=points$bg,alpha=0.6,lwd=points$lwd,cex=points$cex)
    ltext(points$x,points$y,points$lab,cex=points$cex.lab,col=points$col.lab) #,adj=c(0.5,-1))
  }
  
  ## Add latitude/longitude lines
  if (latlon) {
  
    latticks <- c(60,70,80); lonticks <- c(-75,-60,-45,-30,-15)
    panel.contourplot(xx, yy, pdat$lat,subscripts=ii0,lty=2,col="grey40",lwd=0.5,
                      at=latticks,
                      region = FALSE, contour = TRUE, labels = FALSE)
    panel.contourplot(xx, yy, pdat$lon,subscripts=ii0,lty=2,col="grey40",lwd=0.5,
                      at=lonticks,
                      region = FALSE, contour = TRUE, labels = FALSE)
    
    ## Latitude ticks
    xs <- numeric(length(latticks))+min(pdat$xx)
    ys <- c(-3305,-2120,-835)
    yslab <- c("60°N","70°","80°")
    ysrt  <- c(-12,-15,-42)
    for (i in 1:length(ys)) ltext(xs[i],ys[i],yslab[i],srt=ysrt[i],adj=c(-0.1,0),col=1,cex=1.1)
    
    ## Longitude ticks
    xs <- c(-455,-240,-55,110,320)
    ys <- numeric(length(lonticks))+max(pdat$yy)-60
    yslab <- c("-75°","-60°","-45°","-30°","-15°E")
    for (i in 1:length(ys)) ltext(xs[i],ys[i],yslab[i],adj=c(0.5,0),col=1,cex=1.1)
  
  }
  
  # Panel label
  if (!is.na(plabel[1])) {
    x1 <- min(pdat$xx)+130; y1 <- max(pdat$yy)-100
    ltext(x1,y1,plabel[pnow],cex=2.5,col="grey40") #,vfont=c("serif","bold"))
  }
  
  # And title
  if (!is.na(ptext[1])) ltext(300,-2750,ptext[ min(pnow,length(ptext)) ],cex=1.5,pos=1)

}

get.contours <- function(nm,var=NULL,zrange=NULL,n=15,alpha=100,darker=0,
                col=NULL,bias=NULL,at=NULL,at.mark=NULL,rev=FALSE)
{
  
  bias1 <- 1
  
  if (!is.null(var)) {
    zrange1 <- range(abs(var),na.rm=TRUE); zrange1 <- c(0,zrange1[2])
  }
  
  if ( nm == "tt" ) {
    
    at1       <- c(-40,-30,-20,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10)
    at.mark1  <- as.character(at1)
    bias1     <- sum(at1>0) / sum(at1<0)
    col1   <- c("magenta","blue","cyan","white","yellow","red","darkred") 
  
  } else if ( nm %in% c("pp","snow") ) {
    
    if (is.null(var)) zrange1 <- c(0,2500)

    at1      <- pretty(zrange1,n)
    at.mark1 <- as.character(at1)
    col1     <- c("magenta","blue","cyan","green","yellow","red","darkred")
    #col1     <- c("white","yellow","green","darkgreen","cyan","blue","magenta","darkviolet")  
    #col1     <- c("white","yellow","green","cyan","blue","magenta")
    
    cat(" -used pp/snow contours.\n")
    
  } else if ( nm == "smb" ) {
    
    at1      <- c(-6000,-4000,-2000,-1000,-800,-600,-400,-200,-100,0,100,200,400,600,800,1000,2000,4000)
    at.mark1 <- as.character(at1)
    bias1    <- 8/9
    col1     <- c("magenta","blue","cyan","white","yellow","red","darkred")
    #col1     <- c("darkred","red","yellow","white","cyan","blue","magenta")
    
  } else if ( nm %in% c("melt","runoff") ) {
    
    at1      <- c(0,100,200,400,600,800,1000,2000)
    at.mark1 <- as.character(at1)
    col1     <- c("magenta","blue","cyan","white","yellow","red","darkred")
    #col1     <- c("darkred","red","yellow","white","cyan","blue","magenta")
    
  } else if ( nm == "dsmb" ) {
    
    at1      <- c(-6000,seq(-1000,1000,by=100),6000)
    at.mark1 <- as.character(at1)
    at.mark1[1] <- at.mark1[length(at.mark1)] <- ""
    bias1    <- 8/9
    col1     <- c("magenta","blue","cyan","white","yellow","red","darkred")
    #col1     <- c("darkred","red","yellow","white","cyan","blue","magenta")
    
  } else if ( nm %in% c("zs") ) {
    
    at1      <- seq(from=0,to=3300,by=300)
    at.mark1 <- as.character(at1)
    col1     <- c("grey80","white")
    #col1     <- c("darkslategray","white")
    
  } else if ( nm %in% c("land") ) {
    
    at1      <- seq(from=-100,to=3300,by=100)
    at.mark1 <- as.character(at1)
    #at.mark1[which(!at.mark1 %in% c("0","500","1000","1500","2000","2500","3000","3500"))] <- " "
    bias1    <- 5
    col1     <- c("tan","brown")
  
  } else if ( nm %in% c("ocean") ) {
    
    at1      <- seq(from=-4500,to=100,by=100)
    at.mark1 <- as.character(at1)
    #at.mark1[which(!at.mark1 %in% c("0","500","1000","1500","2000","2500","3000","3500"))] <- " "
    col1     <- c("paleturquoise","paleturquoise")
    #col1     <- c("white","white")
  
  } else if ( nm %in% c("dttp1") ) {  # Generic (positive or negative only)

    zrange1  <- zrange
    at1      <- pretty(zrange1,n)
    at.mark1 <- as.character(at1)
    col1     <- c("magenta","cyan","white","yellow","orange","red","darkred")
    
    cat("dttp colors for ",nm,":",zrange1,"\n")
    
  } else if ( nm %in% c("pos","neg") ) {  # Generic (positive or negative only)

    zrange1  <- zrange
    at1      <- pretty(zrange1,n)
    at.mark1 <- as.character(at1)
    col1     <- c("magenta","blue","cyan","yellow","red","darkred")
    
    cat("Generic colors for ",nm,":",zrange1,"\n")
    
  } else {  # Generic (negative and positive)
    
    #zrange1  <- max(abs(zrange))
    #zrange1  <- c(-zrange1,zrange1)
    zrange1  <- zrange
    at1      <- pretty(zrange1,n)
    at.mark1 <- as.character(at1)
    bias1    <- sum(at1>0) / max(1,sum(at1<0))
    col1     <- c("magenta","blue","cyan","lightgreen","white","yellow","orange","red","darkred")
    
    cat("Generic colors for ",nm,":",zrange1,"\n")
  }
  
  # Assign newly generated values if needed
  if (is.null(col))       col     <- col1
  if (is.null(bias))      bias    <- bias1
  if (is.null(at))        at      <- at1
  if (is.null(at.mark))   at.mark <- at.mark1
  
  if (rev) col <- rev(col)
  
  # Now generate actual colors via colorwheel
  # Shift colors darker or lighter, and add alpha layer (transparency)
  colorwheel <- colorRampPalette(col,bias=bias)
  cols <- colorwheel(length(at)-1)
  cols <- darker(cols,darker) 
  cols <- alpha(cols,alpha)
  
  return(list(at=at,at.mark=at.mark,cols=cols,colorwheel=colorwheel,
              contour.plot=FALSE,contour.lwd=0.5,contour.lty=1,contour.col="white"))
}

fill_mask <- function(mask)
{ # Fill an ice sheet mask so there are no isolated grid types
  
  nx <- dim(mask)[1]
  ny <- dim(mask)[2]
  
  fill <- 0
  dn   <- 1
  
  for ( i in 3:(nx-2) ) {
    for ( j in 3:(ny-2) ) {
      
      # For now only fill land points that should be ice
      if ( mask[i,j] == 1 ) {
        
        neighbs <- mask[(i-dn):(i+dn),(j-dn):(j+dn)]
        tot <- base::sum(neighbs)-1
        
        # If all neighbors == 0, then it is isolated by ice
        if ( tot < 4 ) {
          mask[i,j] <- 0
          fill <- fill+1
        }

      }
    }
  }
  
  cat("Filled mask",fill,"times.\n")
  
  return(mask)
}

## Add latitude and longitude lines and labels
add_latlon <- function(grid,x=grid$xx[,1],y=grid$yy[1,],col="grey30",lty=2,lwd=0.5,cex=0.9,shift.lat=0,shift.lon=0)
{
  lats <- grid$lat
  lons <- grid$lon
  
  # Make sure lat/lon is range [-180:180]
  if (max(lons) > 180 ) {
    ii = which(lons>180)
    lons[ii] = lons[ii] - 360
  } 

  # Lat/lons
  latticks <- c(60,70,80)
  lonticks <- c(-75,-60,-45,-30,-15)
  contour(x=x,y=y,z=lats,add=TRUE,levels=latticks,drawlabels=FALSE,lty=lty,lwd=lwd,col=darker(1,30))
  contour(x=x,y=y,z=lons,add=TRUE,levels=lonticks,drawlabels=FALSE,lty=lty,lwd=lwd,col=darker(1,30))

  ## Latitude ticks
  xs <- numeric(length(latticks))+min(x)
  ys <- c(-3315,-2125,-845) - 10 + shift.lat
  yslab <- c("60°N","70°","80°")
  ysrt  <- c(-12,-15,-38)
  for (i in 1:length(ys)) text(xs[i],ys[i],yslab[i],srt=ysrt[i],adj=c(-0.1,0),col=col,cex=cex)
  
  ## Longitude ticks
  xs <- c(-455,-240,-55,110,320)
  ys <- numeric(length(lonticks))+max(y)-72 + shift.lon
  yslab <- c("-75°","-60°","-45°","-30°","-15°E")
  for (i in 1:length(ys)) text(xs[i],ys[i],yslab[i],adj=c(0.5,0),col=col,cex=cex)

}

  
# My own image/contour plotting function
my.image <- function(dat,nm="zs",ii=c(1:length(dat[[nm]])),col=NULL,breaks=NULL,mask=TRUE)
{
  
  x <- dat$xx[,1]; y <- dat$yy[1,]
  
  ## Plot the elevation: filled contour
  var <- dat[[nm]]*NA
  var[ii] <- dat[[nm]][ii]
  
  image(x=x,y=y,z=var,axes=FALSE,col=col,breaks=breaks,xlab="",ylab="")
  
  ## Add land/ice mask
  if (mask) contour(x=x,y=y,z=dat$mask,add=TRUE,drawlabels=FALSE,nlevels=2,lwd=2,col=alpha(1,80))
  contour(x=x,y=y,z=dat$zs,add=TRUE,drawlabels=FALSE,nlevels=10,lwd=0.5,lty=1,col=alpha(1,50))
  
  ## Add lat/lon
  contour(x=x,y=y,z=dat$lat,add=TRUE,drawlabels=FALSE,nlevels=4,lwd=0.5,lty=2)
  contour(x=x,y=y,z=dat$lon,add=TRUE,drawlabels=FALSE,nlevels=4,lwd=0.5,lty=2)
  
  box()
}

my.image3 <- function(dat,nm="zs",suffix=NA,sim=1,k=1,ii=c(1:length(dat[[nm]][sim,,,k])),col=NULL,breaks=NULL,
                      plt.var=TRUE,plt.var.cont=FALSE,plt.ocean=TRUE,plt.land=TRUE,plt.mask=TRUE,plt.latlon=TRUE,lwd=0.5,
                      col.land=NULL,col.latlon=alpha(1,50),cex.latlon=0.9,lwd.mask=2,fill=FALSE)
{
  
  vnm      <- nm
  vnm.mask <- "mask"
  vnm.zs   <- "zs"
  if (!is.na(suffix)) {
    vnm      <- nm
    vnm.mask <- "mask"
    vnm.zs   <- "zs"
  }

  ## Plot the variable: filled contour
  var0 <- dat[[vnm]][sim,,,k]
  var <- var0*NA
  var[ii] <- var0[ii]

  # Get x and y vectors for plotting
  dims = dim(dat[[vnm]][sim,,,k])
  x = seq(0,1,length.out=dims[1])
  y = seq(0,1,length.out=dims[2])
  if ( "xx" %in% names(dat) ) x = dat$xx[,1]
  if ( "yy" %in% names(dat) ) y = dat$yy[1,]
  
  ## Limit variable value to range of contours
  if (!is.null(breaks)) {
    brange <- range(breaks)
    var[var>brange[2]] <- brange[2]
    var[var<brange[1]] <- brange[1]
  }
  
  # Get generic variables
  mask  <- dat[[vnm.mask]][sim,,,k]
  if (fill) mask <- fill_mask(mask)
  
  land  <- dat[[vnm.zs]][sim,,,k]; land[mask != 1]  <- NA
  ocean <- dat[[vnm.zs]][sim,,,k]; ocean[mask != 2] <- NA
  
  zs    = dat[[vnm.zs]][sim,,,k]

#   if (is.null(breaks)) {
#     conts <- get.contours(nm="pos",zrange=range(var,na.rm=TRUE),n=10)    
#     breaks <- conts$at
#     col    <- conts$cols
#   }
  
  # Start the empty plot
  par(xaxs="i",yaxs="i")
  plot(range(x),range(y),type="n",axes=FALSE,xlab="",ylab="")
  
  # Add ocean
  if ( plt.ocean ) {
    cont.ocean <- get.contours(nm="ocean",darker=20)
    ocean[ocean<min(cont.ocean$at)] = min(cont.ocean$at)
    ocean[ocean>max(cont.ocean$at)] = max(cont.ocean$at)
    image(x=x,y=y,z=ocean,add=TRUE,col=cont.ocean$cols,breaks=cont.ocean$at)
  }
  
  # Add land
  if ( plt.land ) {
    cont.land <- get.contours(nm="land",col=col.land)
    land[land<min(cont.land$at)] = min(cont.land$at)
    land[land>max(cont.land$at)] = max(cont.land$at)
    image(x=x,y=y,z=land,add=TRUE,col=cont.land$cols,breaks=cont.land$at)
  }
  
  # Get contours for ice elevation (for consistency
  cont.zs <- get.contours(nm="zs")
  cont.zs$lwd <- rep(lwd,length(cont.zs$at))
  cont.zs$lwd[cont.zs$at %in% c(2400)] <- lwd*2
  
  # Add variable
  if (plt.var) {
    var[var<min(breaks)] = min(breaks)
    var[var>max(breaks)] = max(breaks)
    image(x=x,y=y,z=var,add=TRUE,col=col,breaks=breaks)
    if (plt.var.cont) contour(x=x,y=y,z=var,add=TRUE,col="grey70",lwd=0.8,levels=breaks,drawlabels=FALSE)
  }
  
  ## Add land/ice mask
  if (plt.mask) contour(x=x,y=y,z=mask,add=TRUE,drawlabels=FALSE,nlevels=2,
                        lwd=lwd.mask,col=alpha(1,70))
  
  ## Add land contours
  if (plt.var) {
    contour(x=x,y=y,z=zs,add=TRUE,drawlabels=FALSE,levels=cont.zs$at,
            lwd=cont.zs$lwd,lty=1,col=alpha(1,50))
  }
  
  if (plt.latlon) {
    ## Add lat/lon
    add_latlon(dat,x=x,y=y,col=col.latlon,cex=cex.latlon,shift.lat=30,shift.lon=10)
  }
  
  
  box()
  
  return(list(col=col,breaks=breaks))
}

list_sub <- function(new,default)
{
  if (!is.null(new)) {
    nms1 = names(new)
    for (q in 1:length(nms1)) {
      nm = nms1[q]
      default[[nm]] = new[[nm]]
    }
  }

  return(default)
}

my.image4 <- function(x=NULL,y=NULL,var,zs,mask,lat,lon,ii=c(1:length(var)),col=NULL,breaks=NULL,
                      plt.var=NULL,plt.land=NULL,plt.ocean=NULL,plt.mask=NULL,plt.latlon=NULL )

{ # var : 2D matrix of variable of interest
  # dat : list containing all other background variables
  
  # Get option lists
  plt.var    = list_sub(plt.var,   list(fill=TRUE,cont=FALSE) )
  plt.land   = list_sub(plt.land,  list(fill=TRUE,cont=TRUE,col=alpha(1,50),lwd=0.5) )
  plt.ocean  = list_sub(plt.ocean, list(fill=TRUE,cont=FALSE,col=NULL) )
  plt.mask   = list_sub(plt.mask,  list(fill=TRUE,cont=TRUE,lwd=2,col=alpha(1,70)) )
  plt.latlon = list_sub(plt.latlon,list(cont=TRUE,lwd=0.5,col=alpha(1,50),cex=0.9) )
  
  # Delete unwanted points from var
  var[which(! c(1:length(var)) %in% ii)] = NA 

  # Get x and y vectors for plotting
  dims = dim(var)
  if (is.null(x)) x = seq(0,1,length.out=dims[1])
  if (is.null(y)) y = seq(0,1,length.out=dims[2])

  # Get generic variables
  if (plt.mask$fill) mask <- fill_mask(mask)
  
  land  <- zs; land[mask != 1]  <- NA
  ocean <- zs; ocean[mask != 2] <- NA

#   if (is.null(breaks)) {
#     conts <- get.contours(nm="pos",zrange=range(var,na.rm=TRUE),n=10)    
#     breaks <- conts$at
#     col    <- conts$cols
#   }
  
  # Start the empty plot
  par(xaxs="i",yaxs="i")
  plot(range(x),range(y),type="n",axes=FALSE,xlab="",ylab="")
  
  # Add ocean
  if ( plt.ocean$fill ) {
    cont.ocean <- get.contours(nm="ocean",darker=20)
    ocean[ocean<min(cont.ocean$at)] = min(cont.ocean$at)
    ocean[ocean>max(cont.ocean$at)] = max(cont.ocean$at)
    image(x=x,y=y,z=ocean,add=TRUE,col=cont.ocean$cols,breaks=cont.ocean$at)
  }
  
  # Add land
  if ( plt.land$fill ) {
    cont.land <- get.contours(nm="land",col=plt.land$col)
    land[land<min(cont.land$at)] = min(cont.land$at)
    land[land>max(cont.land$at)] = max(cont.land$at)
    image(x=x,y=y,z=land,add=TRUE,col=cont.land$cols,breaks=cont.land$at)
  }
  
  # Get contours for ice elevation (for consistency
  cont.zs <- get.contours(nm="zs")
  cont.zs$lwd <- rep(plt.land$lwd,length(cont.zs$at))
  cont.zs$lwd[cont.zs$at %in% c(2400)] <- plt.land$lwd*2
  
  # Add variable
  var[var<min(breaks)] = min(breaks)
  var[var>max(breaks)] = max(breaks)
  if (plt.var$fill) image(x=x,y=y,z=var,add=TRUE,col=col,breaks=breaks)
  if (plt.var$cont) contour(x=x,y=y,z=var,add=TRUE,col="grey70",lwd=0.8,levels=breaks,drawlabels=FALSE)
  
  ## Add land/ice mask
  if (plt.mask$cont) contour(x=x,y=y,z=mask,add=TRUE,drawlabels=FALSE,nlevels=2,
                             lwd=plt.mask$lwd,col=plt.mask$col)
  
  ## Add land contours
  if (plt.land$cont) contour(x=x,y=y,z=zs,add=TRUE,drawlabels=FALSE,levels=cont.zs$at,
                             lwd=cont.zs$lwd,lty=1,col=alpha(1,50))
  
  ## Add lat/lon contours
  if (plt.latlon$cont) add_latlon(list(lat=lat,lon=lon),x=x,y=y,col=plt.latlon$col,cex=plt.latlon$cex,shift.lat=30,shift.lon=10)
  
  box() 

  return(list(col=col,breaks=breaks))
}

plot.series <- function(x,y,xlim=NULL,ylim=NULL,col=NULL,good=NULL,
                        axis=c(1,2),lwd=2,lty=1,alpha=20,title=NULL)
{ # x: x variable (vector)
  # y: y variables (array, n X nx)

  if (is.null(xlim)) xlim = range(x,na.rm=T)
  if (is.null(ylim)) {
    kk = which(x >= xlim[1] & x <= xlim[2])
    ylim = range(y[,kk],na.rm=T)

  }
  
  plot(xlim,ylim,type="n",axes=F,ann=F)
  grid()
  axis(1)
  axis(2)
  
  sims = c(1:dim(y)[1])
  if (is.null(good)) good = sims
  
  # First plot bad sims with transparency
  qq = which(! sims %in% good)
  for (q in qq) lines(x,y[q,],lwd=lwd,lty=lty,col=alpha(col[q],alpha))
  
  # Then plot good sims
  qq = which(sims %in% good)
  for (q in qq) lines(x,y[q,],lwd=lwd,lty=lty,col=col[q])
  
  # Plot label
  if (!is.null(title)) text(xlim[1],ylim[2]-diff(ylim)*0.05,pos=4,title,cex=1)
  
  box()
}

load_core_info <- function(filename)
{

  info = read.table(filename,header=TRUE)

  core_info = list()
  for (q in 1:dim(info)[1]) {
    core_info[[info$name[q]]] = info[q,2:dim(info)[2]]
  }

  core_info$table = info 
  
  return(core_info)
}

## Function to plot borehole temperature profiles
plot.core <- function(core,cnm="grip",nm="temp",xlab="Temperature (°C)",ylab="Depth (km)",
                           col=1,bad=NA,ii=NA,lwd=1,convert.y=1e-3,
                           ylim=NULL,xlim=NULL,title=""  )
{
    
  # Define the names (to include the core name)
  ynm <- paste("grip","zs",sep=".")
  xnm <- paste(cnm,nm,sep=".")
  
  # Determine number of runs to plot (and indices of runs)
  nruns <- dim(core[[xnm]])[2]
  if (is.na(ii[1])) ii <- c(1:nruns)
  if (is.na(bad[1])) bad <- numeric(nruns)
  if (length(col)==1) col <- rep(col,nruns)
  if (length(lwd)==1) lwd <- rep(lwd,nruns)
  
  depth <- core[[ynm]]*convert.y
  if (is.null(ylim)) ylim <- range(depth,na.rm=TRUE)
  if (is.null(xlim)) xlim <- range(core[[xnm]],na.rm=TRUE)
  
  # Make the initial (empty) plot
  plot(xlim,ylim,type="n",xlim=xlim,ylim=ylim,xlab="",ylab="")
  #xx <- c(rep(-50,2),rep(10,2)); yy <- c(-4010,10,10,-4010)
  #g1 <- "grey95"; polygon(xx,yy,border=g1,col=g1)
  
  #xx <- c(rep(vline[1],2),rep(vline[2],2)); yy <- c(hline,hline[2:1])
  #polygon(xx,yy,border=NA,col=0)
  grid(col="grey80")
  
  # Loop over runs, but plot grey (bad) runs first
  for ( i in ii ) {
    if ( bad[i] > 0 ) lines(core[[xnm]][,i],depth,lwd=lwd[i],col=col[i])
  }
  
  for ( i in ii ) {
    if ( bad[i] == 0 ) lines(core[[xnm]][,i],depth,lwd=lwd[i],col=col[i])
  }

  #polygon(xx,yy,border=1,lwd=2)
  #abline(h=hline,v=vline,col=2,lwd=2)
  mtext(side=1,line=2,xlab,cex=1)
  mtext(side=2,line=2.5,ylab,cex=1,las=0)
  
  title(main=title)
  
  box()
}

plot2d <- function(mask,zs,mar=NA,plt=NA,title="",cex.title=1,asp=0.7,title.loc=c(0.9,0.15),lwd=1,
                   col=c("white","grey40","grey60"),add=TRUE,lab=TRUE,labcex=0.8,interp=FALSE,shade=FALSE)
{
  def.par <- par(no.readonly=TRUE)
  if ( is.na(mar[1]) ) mar <- def.par$mar
  if ( is.na(plt[1]) ) plt <- def.par$plt
  
  #par(new=add,mar=mar)
  par(new=add,plt=plt)
  
  if (interp != FALSE) {
    # tmp <- list(mask=mask,zs=zs)
    # new <- interp.clim(tmp,nms=c("mask","zs"),factor=2)
    # zs   <- new$zs
    # mask <- new$mask

    zs   = grid.interp(zs,factor=interp)
    mask = grid.interp(mask,factor=interp,mask=TRUE)

  }

  zs[mask==2] <- NA
  image(mask,axes=FALSE,col=col,asp=1/asp)

  #contour(mask,add=TRUE,levels=c(2),drawlabels=FALSE,col=darker(col[2],-5),lwd=lwd)

  contour(zs*1e-3,add=TRUE,levels=seq(from=0.0,to=3.5,by=0.5),lwd=lwd,
          drawlabels=FALSE,col="grey40" )
  
  contour(zs*1e-3,add=TRUE,levels=c(1.5,2.5),lwd=lwd*2.5,
          drawlabels=FALSE,col="grey40")

  #contour(mask,add=TRUE,levels=c(1),drawlabels=FALSE,col=alpha(1,50),lwd=lwd*0.5)
  
  if (shade) add_shade(zs)

  text(title.loc[1],title.loc[2],title,cex=cex.title,pos=2)
  
  if (add) par(new=add,plt=def.par$plt,mar=def.par$mar)  # return to normal
  
}

plot2d.ani <- function(dat, ...) 
{
  cat("Generating frames:",length(dat$time),"..")
  
  for ( i in 1:length(dat$time) ) {
    
    par(mar=c(0,0,0,0))
    plot2d(dat$mask[,,i],dat$zs[,,i],title=paste(dat$time[i],"ka"),add=FALSE )
    
    Sys.sleep(ani.options("interval"))
  }
  
  cat("done.\n")
  
}
