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
