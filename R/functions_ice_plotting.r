
#' @export
plot_icesheet <- function(x,y,z,mask,zb=z*0,add=FALSE,col.cont=alpha("darkcyan",50),
                          xlim=range(x),ylim=range(y))
{
    zs = z
    zs[zs<0] = NA 

    zs_ocean = zb
    zs_land  = z
    zs_ice   = z 
    zs_ocean[mask!=2] = NA  
    zs_land[mask!=1]  = NA  
    zs_ice[mask!=0]   = NA 

    breaks_ocean = pretty(range(zs_ocean,na.rm=TRUE),50)
    cols_ocean   = colorRampPalette(c("grey60","grey95"),bias=0.3)(length(breaks_ocean)-1)

    breaks_land = seq(0,3000,by=50)
    # cols_land   = colorRampPalette(c("wheat3","wheat2","wheat1"))(length(breaks_land)-1)
    cols_land   = colorRampPalette(c("grey95","grey30","black"),bias=1)(length(breaks_land)-1)
    # cols_land   = colorRampPalette(c("darkolivegreen","olivedrab","wheat3","wheat2","wheat4","brown","white"))(length(breaks_land)-1)
    col.cont1   = "grey30"

    breaks_ice = pretty(range(zs_ice,na.rm=TRUE),50)
    cols_ice   = colorRampPalette("white")(length(breaks_ice)-1)

    if (!add) {
        par(plt=c(0,1,0.02,1),xaxs="i",yaxs="i")
        plot(xlim,ylim,type='n',ann=FALSE,axes=FALSE,asp=0.9)
    }

    # image(topo$x,topo$y,zs_ocean,add=TRUE,breaks=breaks_ocean,col=cols_ocean)
    image(x,y,zs_land,add=TRUE,breaks=breaks_land,col=cols_land)
    image(x,y,zs_ice,add=TRUE,breaks=breaks_ice,col=cols_ice)

    add_shade(x,y,z=zs,max=10)

    # contour(topo$x,topo$y,zs_ocean,add=TRUE,col=alpha("darkblue",20),lwd=0.5,drawlabels=F,
    #           levels=c(-5000,-3000,-2000,-1000,-500))
      
    contour(x,y,zs_ice,add=TRUE,levels=seq(0,3500,by=250),drawlabels=FALSE,col=col.cont,lwd=0.5)
    contour(x,y,zs_ice,add=TRUE,levels=c(1000,2000),drawlabels=FALSE,col=col.cont,lwd=1.6)
    contour(x,y,zs_land,add=TRUE,levels=seq(0,2500,by=50),drawlabels=FALSE,col=col.cont1,lwd=0.5)
    contour(x,y,zs_land,add=TRUE,levels=c(1000,2000),drawlabels=FALSE,col=col.cont1,lwd=1.2)
    contour(x,y,mask,  add=TRUE,levels=c(1),drawlabels=FALSE,col=alpha("grey50",10),lwd=5)

}

#' @export
points_cores <- function(corep,cex=corep$cex)
{
    i0 = which(corep$name == "neem")
    i1 = which(corep$name == "neemup")
    if (length(c(i0,i1))==2) 
        segments(corep$x[i0],corep$y[i0],corep$x[i1],corep$y[i1],col="black",lwd=1.5)
    points(corep$x,corep$y,pch=corep$pch,cex=cex,col=corep$col,lwd=corep$lwd,bg=corep$bg)
}
