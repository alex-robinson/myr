library(RNetCDF)

#' @export
my.read.nc = function(filename)
{
    nc = open.nc(filename)
    dat = read.nc(nc)
    close.nc(nc)
    return(dat)
}
