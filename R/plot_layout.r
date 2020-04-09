
#' @export
myfigure <- function(fldr=".",file="Rplot",date=TRUE,type="pdf",engine="cairo",
                     width=NULL,height=NULL,units="mm",asp=1,pointsize=12,res=300,
                     cex=1,cex.lab=1,cex.axis=1,bg="white",onefile=TRUE)
{
    # Some system settings
    host  = system("hostname",intern=TRUE)
    os    = system("uname",intern=TRUE)
    today = format(Sys.time(),"%Y-%m-%d")

    # Make filename
    file = paste(file,".",type,sep="")
    if (date == TRUE) file = paste(today,"_",file,sep="")
    file = file.path(fldr,file)

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

#' @export
mylegend_internal <- function(breaks,col,units="mm",x=c(0,1),y=c(0,1),at=breaks,labels=NULL,
                     xlab="",ylab="",xlim=NULL,ylim=NULL,zlim=range(breaks),
                     cex=1,cex.lab=1,new=FALSE,extend=FALSE,vertical=TRUE,line=1.8,
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

    if (extend != 0) {
      # Add triangles to ends of legend to indicate extended range
        if (vertical) {
            par(xpd=TRUE)
            if (extend %in% c(2,-1)) polygon(c(0,1,0.5),c(0,0,0-0.05),border=col.axis,col=col[1])
            if (extend %in% c(2, 1)) polygon(c(0,1,0.5),c(1,1,1+0.05),border=col.axis,col=col[n-1])
            par(xpd=NA)
        } else {
            par(xpd=TRUE)
            if (extend %in% c(2,-1)) polygon(c(0,0,0-0.05),c(0,1,0.5),border=col.axis,col=col[1])
            if (extend %in% c(2, 1)) polygon(c(1,1,1+0.05),c(0,1,0.5),border=col.axis,col=col[n-1])
            par(xpd=NA)
        }
    }

    par(new=TRUE,xpd=NA,xaxs="i",yaxs="i",...)
    plot(xlim,ylim,type="n",axes=F,ann=F,cex=cex)
    axis(ax,at=at,labels=labels,mgp=mgp,tcl=-0.1,col=col.axis,col.axis=col.axis,cex.axis=cex)
    box(col="grey10")

    mtext(side=1,line=line,xlab,cex=cex.lab)
    mtext(side=2,line=line,ylab,cex=cex.lab)

    par(xpd=FALSE)
}

#' @export
mylegend = function(breaks,col,...,labels=paste(breaks),extend=0,evenspacing=FALSE)
{
    
    if (evenspacing) {
        breaks0 = seq(from=0,to=1,length.out=length(breaks))
        mylegend_internal(breaks=breaks0,col=col,at=breaks0,labels=labels,extend=extend,...)
    } else {
        mylegend_internal(breaks,col,labels=labels,extend=extend,...)
    }
}
#' @export
my.par  <- function(mar=c(3.2,3.3,1,1),xaxs="i",yaxs="i",tcl=0.4,mgp=c(2.5,0.3,0),las=1,...)
{
  par(...,mar=mar,tcl=tcl,mgp=mgp,las=las,xaxs=xaxs,yaxs=yaxs)
}

#' @export
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

#' @export
myplot <- function(...,axes=c(1,2,3,4),box=TRUE,grid=TRUE,
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

#' @export
plotblank <- function(mar=NA)   
{ ## plot an empty space, in which text can be added
  #def.par <- par(no.readonly=TRUE)
  #if (is.na(mar[1])) mar <- c(1,1,1,1)
  #par(mar=mar)
  plot(c(0,1),c(0,1),type="n",axes=FALSE, xlab="", ylab="")
  
  #par(def.par$mar) # return to normal
}
