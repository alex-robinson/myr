
# Depends on filters...
# source("myr/filters.r")

#' @export
corr_oceanic = function(time,d18O_ice,time_ocean=time,d18O_ocean)
{   # Following Jouzel et al. (2003) and Kindler et al. (2014)
    # Corr_ocean = del[dD_ocean] x (1 + dD_ice)/(1 + del[dD_ocean])

    k0  = which(time_ocean==0)
    tmp = d18O_ocean - d18O_ocean[k0]
    del_ocean = approx(time_ocean,d18O_ocean,xout=time)$y 

    k0 = which(time==0)
    del_ice = d18O_ice-d18O_ice[k0]

    corr = d18O_ice - del_ocean * (1+d18O_ice) / (1+del_ocean)
    # corr = d18O_ice - del_ocean * (1+del_ice) / (1+del_ocean)

    return(corr)
}

#' @export
convert_dT   = function(time,d18O,dTlgm,lgm=c(-22,-20))
{   # Convert d18O to temperature anomaly
    
    # new = data.frame(time=xout,d18O=approx(time,d18O,xout=xout,rule=2)$y)

    k0    = which.min(abs(time-0)) 
    kkLGM = which(time >= lgm[1] & time <= lgm[2])

    dT = (d18O-d18O[k0])/(d18O[k0]-mean(d18O[kkLGM],na.rm=TRUE)) * abs(dTlgm)

    return(dT)
}

#' @export
convert_T_kindler = function(time,d18O,alpha,beta,T_pd=241.6,d18O_pd=35.1)
{   # Convert NGRIP d18O to temperature anomaly following fit parameters for 
    # glacial time period
    #
    T = (d18O + d18O_pd)/alpha + T_pd + beta 

    return(T)
}

#' @export
convert_T_jouzel = function(time,dD,alpha,beta)
{   # Convert EPICA Dome C dD to T
    # Jouzel et al. (2007): δD = 6.2 ‰/°C*T+ 5.5‰  == T = (dD-5.5[pmil])1/6.2[pmil/degC]
    #
    T = (dD - beta)/alpha

    return(T)
}


## Function to generate subsurface ocean temperature index
## When grip is in stadial (state 0), subsurface is warming
## When grip is in interstadial (state 1), subsurface is cooling
#' @export
relaxer <- function(time,S,Tmin=0,Tmax=1,taumin=0.1,taumax=1)
{ # 
  dt = diff(time[1:2])
  n  = length(time)
  
  # Initialize vector of temperatures with 0
  T = numeric(n)
  T00 = tau = dTdt = numeric(n)

  for (k in 2:n) {
    
    if (!is.na(S[k])) {

        # Determine constants based on state (cooling or warming)
        if ( S[k] == 0 ) {
          T00[k]  = Tmax
          tau[k] = taumax
        } else {
          T00[k]  = Tmin
          tau[k] = taumin
        }

        dTdt[k] = (T00[k]-T[k-1]) / tau[k]

        T[k] = T[k-1] + dTdt[k]*dt

    }

  }
  
  # Make sure zeros are zeros
  T[ T < 1e-3 ] = 0

  return(T)
}

#' @export
calc_stadials <- function(time,d18O,sdfac=1.2,time0=min(time),pc=c(19,0.1),L=2)
{
    dat = data.frame(time=time,d18O=d18O)

    # Get the millennial signal for detrending
    if (is.null(pc)) {
        dat$d18Omil = dat$d18O 
    } else {
        ii = which(!is.na(dat$d18O))
        dat$d18Omil =dat$d18O*NA
        dat$d18Omil[ii] = myr::freqfilter(x=dat$time[ii],y=dat$d18O[ii],pc=pc,bp=TRUE)
    }

    ii = which(is.na(dat$d18Omil))
    dat$d18Omil[ii] = approx(dat$time,dat$d18Omil,xout=dat$time[ii],rule=2)$y

    # # Smooth the detrended signal
    # # dat$d18Omil.sm = ssatrend(dat$d18Omil,L=20)   # L=20 for dt=50 => 2 ka smoothing window
    # dt = diff(dat$time[1:2])
    # nt = length(dat$time)
    # sm = loess(dat$d18Omil ~ dat$time,span=(2/dt)/nt)
    # dat$d18Omil.sm = predict(sm)
    # tmp0 = dat$d18Omil.sm 

    # Now take 1st derivative
    tmp0 = dat$d18Omil
    dtmp = c(0,diff(tmp0)/diff(dat$time))
    dtmp[1] = dtmp[2] 
    
    # Smooth the 1st derivative
    dt = diff(dat$time[1:2])
    nt = length(dat$time)
    sm = loess(dtmp ~ dat$time,span=(L/dt)/nt)
    dtmp = predict(sm)

    # Get sd and threshold value
    std  = sd(dtmp,na.rm=T)
    threshhi =  std*sdfac
    threshlo = -std*sdfac

    # Generate climate state time series
    dat$state = numeric(length(dat$time))+1      # Starting with interstadial
    new = 0 
    kk = which(dat$time >= time0)
    for (k in kk) {
        if ( is.na(dtmp[k]) ) {
            new = NA 
        } else if ( dtmp[k] <  threshlo ) {
            new = 0
        } else if ( dtmp[k] >  threshhi ) {
            new = 1
        }
        dat$state[k] = new
    }

    # Save derivative
    dat$Dd18Omil = dtmp

    # Generate subsurface ocean temp index
    dat$subsurface = relaxer(dat$time,dat$state,taumin=0.1,taumax=1.0)
  
    kk = which(dat$time < time0)
    dat$state[kk] = NA 
    dat$subsurface[kk] = NA

    return(dat)
}
  
