# For countries and maps
# require(maptools)
data(wrld_simpl)
# require(gpclib)
# require(rgeos)
gpclibPermit()

#' @export
wrld180 = wrld_simpl
wrld360 = nowrapRecenter(wrld_simpl,avoidGEOS=TRUE)

#' @export
world_countries <- function() return(as.vector(wrld_simpl$NAME))

#' @export
mask_by_region <- function(lon,lat,region=NULL)
{   # http://gis.stackexchange.com/questions/75033/given-a-lat-and-lon-identify-if-the-point-is-over-land-or-ocean

    if (is.null(region)) region = expand.grid(lon=lon,lat=lat)

    #region = rbind(region,region[1,])

    # Expand the grid array dimensions (lon,lat) to a set of data points
    points <- expand.grid(lon=lon,lat=lat)

    # Determine which points in the grid are inside of the regional polygon
    ii <- which(sp::point.in.polygon(points[,1],points[,2],region[,1],region[,2]) > 0)

    # Make a mask of points inside the region
    mask = array(0,dim=c(length(lon),length(lat)))
    mask[ii] = 1

    return(mask)
}

#' @export
mask_by_country <- function(lon,lat,countries="United States")
{   # http://gis.stackexchange.com/questions/75033/given-a-lat-and-lon-identify-if-the-point-is-over-land-or-ocean

    lon[lon>=180] = lon[lon>=180]-360 

    ## Create a SpatialPoints object
    points <- expand.grid(lon, lat)  # Note that I reversed OP's ordering of lat/long
    pts <- sp::SpatialPoints(points, proj4string=sp::CRS(sp::proj4string(wrld_simpl)))

    ## Find which points fall over land
    qq = wrld_simpl$NAME %in% countries 

    ii <- !is.na(sp::over(pts, wrld_simpl[qq,])$FIPS)

    mask = array(0,dim=c(length(lon),length(lat)))
    mask[ii] = 1

    return(mask)
}

#' @export
my.in.land.grid <- function(lon,lat)
{ 
  countries = world_countries()
  mask = mask_by_country(lon,lat,countries)

  return(mask)
}


## Regional indices on Earth
#' @export
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
#' @export
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

#' @export
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

#' @export
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

#' @export
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

#' @export
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

#' @export
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

#' @export
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

#' @export
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

