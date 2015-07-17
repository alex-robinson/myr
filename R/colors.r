
# Default jet colors
#' @export
jet.colors = c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F",
                 "yellow", "#FF7F00", "red", "#7F0000")

#' @export
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

#' @export
alpha <- function(col,percent=50)
{
  c <- rep("#000000",length(col))
  for (i in 1:length(col)) {
    cc <- col2rgb(col[i])/255
    c[i] <- rgb(t(cc),alpha=percent/100)
  }
  
  return(c)
}

#' @export
colmix <- function(col,col0="grey80",percent=50) {
    newcol = col
    for (q in 1:length(col)) newcol[q] = colorRampPalette(c(col0,col[q]))(100)[percent]
    return(newcol)
}

#' @export
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

#' @export
get_col = function(x,col=c("blue","red"),n=20,xlim=range(x,na.rm=TRUE))
{
  nx    = length(x)
  dxlim = diff(xlim)

  breaks = pretty(xlim,n)
  db = diff(breaks)[1]
  nb = length(breaks)

  palette = colorRampPalette(col)(nb-1)
  cols = rep("black",nx)
  for ( i in 1:nx ) {
    j = which.min( abs(x[i] - (breaks[1:(nb-1)]+db/2)) )
    if ( length(j)==1 ) cols[i] = palette[j]
  }
  
  return(list(breaks=breaks,palette=palette,col=cols))
}
