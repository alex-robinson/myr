


calc_core_meta = function(core,time_range,conv_time=1)
{
    meta = list()
    nsim = length(core$sim)

    # Save dimension variables
    meta$time   = core$time*conv_time
    meta$age_ka = -meta$time
    meta$depth  = seq(-3800,0,by=50)
    meta$sim    = core$sim

    # Get present-day depth profile 
    temp_p = core$temp_p 
    age_p  = core$age_p 

    # Interpolate to new age scale 
    
    meta$temp_p = array(NA,dim=c(nsim,length(meta$depth)))
    meta$age_p  = array(NA,dim=c(nsim,length(meta$age_p)))
    
    k0 = which.min(abs(meta$time-0.0))
    for (i in 1:nsim) {
        depth = core$zc[i,] - max(core$zc[i,]) 

        if (core$H[i,k0] > 0) {
            meta$temp_p[i,] = approx(depth,temp_p[i,],xout=meta$depth)$y 
            meta$age_p[i,]  = approx(depth,age_p[i,], xout=meta$depth)$y 
        }
    }

    # Which variables need anomalies 
    nms = c("zs","H","tt","tjja","ttp","pp","snow")

    # Get present-day values and anomaly time series
    k0 = which.min(abs(meta$time-0.0))
    meta$today = list(time=meta$time[k0])
    for (q in 1:length(nms)) {
        nm = nms[q]
        meta[[nm]]       = core[[nm]]
        meta$today[[nm]] = core[[nm]][,k0]

        dnm = paste0("d",nm)
        meta[[dnm]] = meta[[nm]] 
        for (i in 1:nsim) meta[[dnm]][i,] = meta[[nm]][i,] - meta$today[[nm]][i]

    }

    # Figure out indices of time range of interest 
    kk = which(meta$time >= time_range[1] & meta$time <= time_range[2])

    # Get minimum levels and times for each variable 
    meta$min = list()
    for (q in 1:length(nms)) {

        nm    = nms[q]
        dnm   = paste0("d",nm)
        nm_tm = paste0(nm,"_time")

        meta$min[[nm]]    = rep(NA,nsim)
        meta$min[[dnm]]   = rep(NA,nsim)
        meta$min[[nm_tm]] = rep(NA,nsim)

        for (i in 1:nsim) {
            k  = kk[which.min(meta[[nm]][i,kk])] 
            meta$min[[nm]][i]    = meta[[nm]][i,k]
            meta$min[[dnm]][i]   = meta[[dnm]][i,k]
            meta$min[[nm_tm]][i] = meta$time[k]
        }
    }

    # Get maximum levels and times for each variable 
    meta$max = list()
    for (q in 1:length(nms)) {

        nm    = nms[q]
        dnm   = paste0("d",nm)
        nm_tm = paste0(nm,"_time")

        meta$max[[nm]]    = rep(NA,nsim)
        meta$max[[dnm]]   = rep(NA,nsim)
        meta$max[[nm_tm]] = rep(NA,nsim)

        for (i in 1:nsim) {
            k  = kk[which.max(meta[[nm]][i,kk])] 
            meta$max[[nm]][i]    = meta[[nm]][i,k]
            meta$max[[dnm]][i]   = meta[[dnm]][i,k]
            meta$max[[nm_tm]][i] = meta$time[k]
        }
    }

    # Get duration of ice-free conditions 
    meta$H0_t0 = rep(NA,nsim)
    meta$H0_t1 = rep(NA,nsim)
    meta$H0_dt = rep(0,nsim)

    for (i in 1:nsim) {
        if (min(meta$H[i,kk])==0) {
            kk1 = kk[which(meta$H[i,kk]==0)]
            meta$H0_t0[i] = meta$time[kk1[1]]
            meta$H0_t1[i] = meta$time[kk1[length(kk1)]]
            meta$H0_dt[i] = meta$H0_t1[i] - meta$H0_t0[i] 
        }
    }

    return(meta)
}

calc_series_meta = function(time,var,time_range)
{
    nsim = dim(var)[1]

    out        = list(time=time,var=var)
    out$age_ka = -out$time
    
    # Also get anomalies relative to present 
    k0 = which.min(abs(time-0.0))
    if (abs(time[k0]) < 0.1) {
        out$dvar = var 
        for (i in 1:nsim) out$dvar[i,] = out$var[i,] - out$var[i,k0]
    }

    meta = data.frame(sim=c(1:nsim))
    meta$min      = NA 
    meta$min_time = NA 
    meta$max      = NA 
    meta$max_time = NA 

    # Figure out indices of time range of interest 
    kk = which(time >= time_range[1] & time <= time_range[2])

    for (i in 1:nsim) {
        k  = kk[which.min(var[i,kk])] 
        meta$min[i]      = var[i,k]
        meta$min_time[i] = time[k] 

        k  = kk[which.max(var[i,kk])] 
        meta$max[i]      = var[i,k]
        meta$max_time[i] = time[k] 
    }

    out$meta = meta 

    return(out)
}

calc_sealevel = function(time,Vtot,dVdt,Vtot0=3.7,z_sle_pd=7.3)
{   
    # Define conversion constant (Gt => sle; Gt/a => slr)
    km32sle = 7.3 / 2.93       # 7.3 [m] / 2.93 [1e6 km^3]
    Gt2slr = ( (1e6/1e3) /3.62e8 ) *1000 # 1e6: giga->kg; 1e3: kg/m3 water; SA: m2 area; 1000: m-> mm
   
    meta = list(time=time,Vtot=Vtot,dVdt=dVdt,Vtot0=Vtot0) 

    # Sea level equivalent, rel. to present, and rate 
    meta$z_sle  = meta$Vtot *km32sle   # m sle 
    meta$z_slr  = z_sle_pd-meta$z_sle  # m sle
    meta$dz_slr = -meta$dVdt*Gt2slr    # (mm/a or m/ka)

    # Different slr calculation!
    meta$z_slrp   = meta$z_slr*NA 
    meta$z_slr2   = meta$z_slr*NA 
    meta$dz_slr2  = meta$z_slr*NA 

    if (length(Vtot0)==1) Vtot0 = rep(Vtot0,dim(meta$Vtot)[1])

    for (i in 1:dim(meta$Vtot)[1]) {
        meta$z_slrp[i,]  = 100*(Vtot0[i]-meta$Vtot[i,])/Vtot0[i]
        meta$z_slr2[i,]  = z_sle_pd*(Vtot0[i]-meta$Vtot[i,])/Vtot0[i]
        meta$dz_slr2[i,] = c(0,diff(meta$z_slr2[i,])/diff(meta$time))
    }

    return(meta)
}

