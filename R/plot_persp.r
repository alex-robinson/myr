
#' @export
mypersp = function(zs,mask,slr=0,time=0,dT=6,col=c("white","wheat2","white"),...)
{
  zs[zs<0] = 0
  
  # Correct mask so only points at sea level are ocean
  #mask[mask==2 & zs > 0] = 1
  
  nx = dim(zs)[1]
  ny = dim(zs)[2]
  x = c(1:nx)
  y = c(1:ny)

  obj.zs   = list(x=x,y=y,z=zs)
  obj.mask = list(x=x,y=y,z=mask)
  
  x1 = seq(1,x[nx],length.out=nx*4)
  y1 = seq(1,y[ny],length.out=ny*4)
  newgrid = list(x=x1,y=y1)

  obj.zs1   = interp.surface.grid(obj=obj.zs,grid.list=newgrid)
  obj.mask1 = interp.surface.grid(obj=obj.mask,grid.list=newgrid)
  obj.mask1$z = round(obj.mask1$z)

  maskcol = drape.color(z=obj.mask1$z,col=col)
  res = persp(x=obj.zs1$x,y=obj.zs1$y,z=obj.zs1$z,theta=20,phi=55,expand=0.004,shade=0.3,
              box=F,d=1,col=maskcol$color.index,border=NA,scale=F,zlim=c(0,1))
  
  pts = trans3d(mean(x),mean(y),3500,pmat=res)
  points(pts,pch=20,col=2)

  par(new=TRUE,mar=rep(0,4))
  plot(c(0,1),c(0,1),type="n",axes=F,ann=F)
  
  ## Add a North arrow
  afac = 0.1
  ax = c(0,0.4)*afac + 0.55
  ay = c(0,0.92)*afac  + 0.1
  arrows(x0=ax[1],y0=ay[1],x1=ax[2],y1=ay[2],angle=30,length=0.2,lwd=2,col="grey30")
  #text(0.58,0.12,"N",srt=-28,cex=1)
  
  # Define a projected N (line by line, 3 line segments total)
  # c(x0,y0,x1,y1)
  nfac = 0.03
  noff = c(0.57,0.1,0.57,0.1)
  n1 = c(-0.2,0,0.2,1)*nfac + noff
  n2 = c(0.2,1,0.2,-0.15)*nfac + noff
  n3 = c(0.2,-0.15,0.6,0.85)*nfac + noff

  segments(x0=n1[1],y0=n1[2],x1=n1[3],y1=n1[4],lwd=2,col="grey30")
  segments(x0=n2[1],y0=n2[2],x1=n2[3],y1=n2[4],lwd=2,col="grey30")
  segments(x0=n3[1],y0=n3[2],x1=n3[3],y1=n3[4],lwd=2,col="grey30")

  # ## Add a title and some additional info
  # text(0.02,0.93,paste("Greenland ice sheet under ",dT,"°C global warming",sep=""),cex=1.2,pos=4)
  # par(family="mono")

  # xstr = format(time*1e3,trim=T,width=4,nsmall=0)
  # ystr = format(round(slr,1),trim=T,width=4,nsmall=1)
  # text(0.03,0.82,paste("  ELAPSED TIME:  ",xstr," years",sep=""),cex=1,col="grey15",pos=4)
  # text(0.03,0.78,paste("SEA LEVEL RISE:  ",ystr," m ",sep=""),cex=1,col="grey15",pos=4)
  
  
}

#' @export
shading <- function(zs,dx=20e3,altitude=pi/6,azimuth=pi/2,scale=1,eps=NULL)
{   # From this python notebook:
    # http://nbviewer.ipython.org/github/ThomasLecocq/geophysique.be/
    #   blob/master/2014-02-25%20Shaded%20Relief%20Map%20in%20Python.ipynb

    dxy      = hgrad(zs,dx=dx)
    slope    = pi/2 - atan(sqrt(dxy$x^2 + dxy$y^2))
    aspect   = atan2(-dxy$x,dxy$y)

    shaded   = sin(altitude) * sin(slope) + cos(altitude) * cos(slope) *
                    cos((azimuth - pi/2) - aspect)

    shaded = shaded-min(shaded,na.rm=TRUE)
    shaded = shaded/max(shaded,na.rm=TRUE)

    # Eliminate points that have no aspect ratio (+/- eps)
    if (!is.null(eps)) {
        shaded[shaded>=(0.5-eps) & shaded <= (0.5+eps)] = NA 
        shaded = shaded*scale
    }
    
    return(shaded)
}
