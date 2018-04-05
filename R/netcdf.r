library(RNetCDF)

#' @export
my.read.nc = function(filename,verbose=TRUE)
{
    nc = open.nc(filename)
    
    if (verbose) {
        info = file.inq.nc(nc)
        cat("ndims=",info$ndims,", nvars=",info$nvars,"\n")

        if (info$unlimdimid > 0) {
            dim_unlim_name = var.inq.nc(nc, info$unlimdimid)$name
            dim_unlim      = var.get.nc(nc,info$unlimdimid)
            cat(dim_unlim_name,": ",dim_unlim,"\n")
        }
    }

    # Now load the data and close the file
    dat = read.nc(nc)
    close.nc(nc)

    if (verbose) {
        cat("vars: ",names(dat),"\n")
    }

    return(dat)
}
