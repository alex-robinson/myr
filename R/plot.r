
#' @export
my.polygon <- function(x,ymin,ymax,...)
{
  if (length(ymin)==1) ymin = rep(ymin,length(x))
  if (length(ymax)==1) ymin = rep(ymax,length(x))
  
  pp = list(x=c(x,rev(x)),y=c(ymin,rev(ymax)))
  polygon(pp,...)
}


## FOR 3D plots with shading ## 

#' @export
col.grad <- function(zs,ltheta=-135,lphi=20)
{ # Generate a scalar to increase or decrease shading
  # -1: darkest color
  #  0: no change
  #  1: lightest color 
  
  # Calculate vector of sunshine (from sun to point)
  # negative, bc points in towards origin...
  torads = pi/180
  dx = sin(lphi*torads)*cos(ltheta*torads)
  dy = sin(lphi*torads)*sin(ltheta*torads)
  dz = -cos(lphi*torads)
  vsun = c(dx,dy,dz)

  # First get surface normals
  vsn = surface.normal2(zs)
  
  # Make into a handy vector for calculating angles
  #vsn = cbind(as.vector(sn$dxx),as.vector(sn$dyy),as.vector(sn$dzz))
  
  # Calculate angle of surface from sunshine for each point
  angles = apply(vsn,FUN=vec.angle,MARGIN=c(1,2),y=vsun)

  # angles = vsn[1,]*0
  # for ( k in 1:dim(vsn)[1] ) {
  #   angles[k] = vec.angle(vsun,vsn[k,])
  # }
  #dim(angles) = dim(zs)
  
  # The larger the angle, the brighter the color should be,
  # 180: directly pointing at surface
  #   0: directly pointing away from surface
  #colshift = (angles/180 - 0.5)
  
  # Get the partially transparent colormask
  # color is black, full transparency or none for brightness
  colshift = 1-(angles/180)
  #colshift = rgb(0,0,0,alpha=colshift)

  return(colshift)
}

#' @export
add_shade <- function(zs,plot=TRUE)
{
  cat("Calculating contour shading...")
  colshift = col.grad(zs)
  colshade=rgb(1,1,1,seq(0,40,length.out=100),maxColorValue=100)
  cat("done.\n")
  
  if (plot) image(colshift,col=colshade,add=T)

  return(list(level=colshift,col=colshade)) 
}

#' @export
test.shading <- function(zs,var=NULL)
{

  shade = add_shade(zs,plot=FALSE)

  elev = get.contours("zs")
  image.plot(a$zs,breaks=elev$at,col=elev$cols)
  image(shade$level,col=shade$col,add=T)
  
  breaks = pretty(range(var,na.rm=T),50)
  cols   = colorRampPalette(jet.colors)(length(breaks)-1)
  cols = alpha(cols,65)
  if (!is.null(var)) image(var,add=T,breaks=breaks,col=cols)
  
  contour(a$zs,add=T,levels=seq(0,3500,by=250),drawlabels=FALSE,col="grey50")
  contour(a$zs,add=T,levels=c(2500),drawlabels=FALSE,col="grey50",lwd=3)

}

###############

### Plotting vectors ###

par.uin <- function()
  # determine scale of inches/userunits in x and y
  # from http://tolstoy.newcastle.edu.au/R/help/01c/2714.html
  # Brian Ripley Tue 20 Nov 2001 - 20:13:52 EST
{
    u <- par("usr")
    p <- par("pin")
    c(p[1]/(u[2] - u[1]), p[2]/(u[4] - u[3]))
}

#' @export
quiver<- function(xx=NULL,yy=NULL,uu,vv,scale=0.13,length=0.008,add=TRUE,col=1,lwd=1,lty=1,thin=1)
# first stab at matlab's quiver in R
# from http://tolstoy.newcastle.edu.au/R/help/01c/2711.html
# Robin Hankin Tue 20 Nov 2001 - 13:10:28 EST
{
    if (thin > 1) {
      # Thin arrays for vector plotting 
      ii = seq(2,dim(xx)[1],by=as.integer(thin))
      jj = seq(2,dim(xx)[2],by=as.integer(thin))
      xx = xx[ii,jj]
      yy = yy[ii,jj]
      uu = uu[ii,jj]
      vv = vv[ii,jj]
    }

    if (is.null(xx)) {
      xx <- col(uu)
      xx = xx/max(xx)
    }
    if (is.null(yy)) {
      yy <- max(row(uu))-row(uu)
      yy = yy/max(yy)
    }    

    speed <- sqrt(uu*uu+vv*vv)
    maxspeed <- max(speed,na.rm=TRUE)

    ux <- uu*scale/maxspeed
    vy <- vv*scale/maxspeed

    #matplot(xpos,ypos,add=add,type="p",cex=0)
    arrows(xx,yy,xx+ux,yy+vy,length=length*min(par.uin()),col=col,lwd=lwd,lty=lty)
}

###############

