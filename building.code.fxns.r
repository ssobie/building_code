##Contains functions to compute Building Code values from the BC
##Building Code and the American ASHRAE Handbook

library(extRemes)
library(ismev)

##------------------------------------------------------------------
##Constants
a.lapse  <<- 0.125      ##km C-1 
cpd      <<- 1004.67    ##J kg-1 K-1 ::Specific Heat Capacity of Dry Air
sat.vape <<- 0.611      ##kPa        ::Saturation vapour pressure at 273 K
grav     <<- 9807.0     ##km s-2     ::Gravitational acceleration
lh.vape  <<- 2.501*10^6 ## J kg-1    ::latent heat of vaporization at 273 K
R.dry    <<- 287.053    ##J K-1 kg-1 ::Gas constant for dry air
R.vapor  <<- 461        ##J K-1 kg-1 ::Gas constant for water vapor
T.zero   <<- 273        ##K          ::Zero in Kelvin
d.lapse  <<- 9.8        ##K km-1     ::dry-adiabatic lapse rate
epsilon  <<- 0.622      ##kgw/kga    ::kg water / kg dry air


##Dew Point Temperature
dew.point.temp <- function(pas,sp.hum) {
     ##Vapor Pressure
     vape.press <- sp.hum * pas / epsilon

     ##Dry bulb temp
     temp.dew <- (1/T.zero - ((R.vapor/lh.vape) * log(vape.press/sat.vape)))^-1
     return(temp.dew - 273)
}

##Wet Bulb Temperature
temp.wet.bulb <- function(temp,temp.dew,pas,sp.hum) {

     ##Mixing Ratio
     mix.ratio <- sp.hum / (1 - sp.hum)

     ##Specific Heat Capacity
     heat.cap <- cpd * (1 + 0.84*mix.ratio) 

     ##Saturation Vapour Pressure
     vp.sat <- sat.vape * exp( (lh.vape/R.vapor) * (1/T.zero - 1/temp))

     ##Saturation Mixing Ratio
     sat.mix.ratio <- (epsilon * vp.sat) / (pas - vp.sat)

     ##Saturated Adiabatic Lapse Rate
     a.lapse.sat <- (grav/heat.cap) * (1 + (sat.mix.ratio * lh.vape)/(R.dry * temp)) / 
                                      (1 + (lh.vape^2 * sat.mix.ratio * epsilon)/(heat.cap * R.dry * temp^2))

     ##Lifting Condensation Level Height
     lift.height <- a.lapse * (temp - temp.dew)

     ##Air temperature at the lifting condensation level
     temp.lcl <- temp - (d.lapse * lift.height)

     ##Wet Bulb Temperature
     temp.WB <- temp.lcl + (a.lapse.sat * lift.height)
     return(temp.WB)
}


##------------------------------------------------------------------
##Driving Rain Wind Pressure
wind.rain.press <- function(uas,vas,pas,tas,huss) {

   R.dry <- 287.053 ## J K-1 kg-1 ::Gas constant for dry air
   tas <- tas+273
   mix.ratio <- huss / (1 - huss)
   tas.virt <- tas * (1 + 0.61*mix.ratio)
   rho.air <- pas / ( R.dry * tas.virt)
   wspd <- sqrt(uas^2 + vas^2)
   
   wind.load <- 0.5 * rho.air * wspd^2
   return(wind.load)
}

##Wind Direction in degrees with North at 0 and East at 90 (clockwise)
wind.dir <- function(uas,vas) {        
   wdir <- (270 - atan2(vas,uas)/pi*180)%%360
   return(wdir)
}


##------------------------------------------------------------------
##Enthalpy

enthalpy <- function(tas,huss) {
   tas <- tas+273      
   cpd   <- 1004.67 ##J kg-1 K-1

   mix.ratio <- huss / (1 - huss)
   heat.cap <- cpd * (1 + 0.84 * mix.ratio)
   enthap <- heat.cap * tas
   return(enthap)  
}


##--------------------------
##Common Parameters

mean.fxn <- function(data,fac){tapply(data,fac,mean,na.rm=T)}
sd.fxn   <- function(data,fac){tapply(data,fac,sd,na.rm=T)}
max.fxn  <- function(data,fac){tapply(data,fac,max,na.rm=T)}
min.fxn  <- function(data,fac){tapply(data,fac,min,na.rm=T)}
sum.fxn  <- function(data,fac){tapply(data,fac,sum,na.rm=T)}
quant.fxn  <- function(data,fac,pctl){tapply(data,fac,quantile,pctl,na.rm=T)}

cold.mon <- function(data,fac){tapply(data,fac,function(x) {length(x[x<0 & !is.na(x)])})} ##Number of days below freezing 
wet.mon  <- function(data,fac){tapply(data,fac,function(x) {length(x[x>1 & !is.na(x)])})} ##Number of days with precip
 
##Degree Day Values
dd <- function(tmean,tbase) {
  g <- tmean - tbase
  days <- sum(g[g>0], na.rm=T)
  return(round(days))
}
cdd.24.fxn  <- function(data,fac){tapply(data,fac, dd, tbase=24)}                            ##Cooling degree days
cdd.18.3.fxn  <- function(data,fac){tapply(data,fac,dd, tbase=18.3)}                           ##Cooling degree days below 18.3
cdd.18.fxn  <- function(data,fac){tapply(data,fac,dd, tbase=18.3)}                           ##Cooling degree days below 18.3
cdd.fxn  <- function(data,fac){tapply(data,fac, dd, tbase=10)}                            ##Cooling degree days
hdd.fxn  <- function(data,fac){tapply(-data,fac,dd, tbase=-10)}                           ##Heating degree days below 10
hdd.18.fxn  <- function(data,fac){tapply(-data,fac,dd, tbase=-18)}                           ##Heating degree days below 18
hdd.18.3.fxn  <- function(data,fac){tapply(-data,fac,dd, tbase=-18.3)}                           ##Heating degree days below 18.3

pctl.fxn <- function(data,fac,pctl){tapply(data,fac,quantile,pctl,na.rm=T)}

conditional.mon.fxn <- function(main,second,fac,pctl) {    
    margin <- tapply(second,fac,quantile,pctl,na.rm=TRUE)                            
    mon.avg <- rep(0,12)

    for (i in 1:length(margin)) {
        ix <- levels(fac)[i]
        ix.subset <- fac %in% ix
        data.subset <- main[ix.subset]
        sec.subset <- second[ix.subset]
        mon.avg[i] <- mean(data.subset[sec.subset > margin[i]],na.rm=T)
    }
    return(mon.avg)
}

##----------------------------------------------------------------------
##Return Periods
calc.return.periods <- function(data,dates,var.name,rperiod) {
  
  yearly.fac <- as.factor(format(dates,'%Y'))
 
  if (sum(is.na(data)) == length(data)) {
    return(NA)
  } else {
    fx <- switch(var.name,
                 tasmax=max,
                 tasmin=min,
                 pr=max,
                 snd=max)
                 
    ts.yearly <- tapply(data,list(yearly.fac),fx,na.rm=T)
    inf.flag <- is.infinite(ts.yearly)
    ts.yearly[inf.flag] <- NA
    
    na.flag <- is.na(ts.yearly)
    if(sum(na.flag)>0) {
      ts.to.fit <- as.vector(ts.yearly[-which(na.flag)])
    } else {
      ts.to.fit <- as.vector(ts.yearly)
    }
    if (var.name=='tasmin') {
      ts.to.fit <- -ts.to.fit
      u.len <- length(unique(ts.to.fit))
      f.len <- length(ts.to.fit)
      if (u.len < (f.len*2/3))
        ts.to.fit <- jitter(ts.to.fit,amount=3)
    }

    if (var.name=='snd') {
      flag <- ts.to.fit == 0
      ts.to.fit[flag] <- jitter(mean(ts.to.fit[!flag]),amount=2)
      ts.to.fit[ts.to.fit<0] <- 1
      print(ts.to.fit)

    }      
        
    ts.fit <- fevd(ts.to.fit,type='GEV')
    ts.old <- gev.fit(ts.to.fit,show=FALSE)

    ts.fit$results$par <- ts.old$mle
    names(ts.fit$results$par) <- c('location','scale','shape')
    ts.rp <- return.level(ts.fit,return.period=as.numeric(rperiod),make.plot=F)

    ##This gives a 3 element vector with lower bound, rp value and upper bound
    rv <- as.numeric(ts.rp)
    print(rv)
    
    if (var.name=='tasmin') {
      rv <- -as.numeric(ts.rp)
    }
    return(rv)
  }
}

get.var.units <- function(var.name) {
  leg.label <- 'degC'
  if (grepl("(pr|rx|r9|RP|rp)", var.name))
    leg.label <- 'mm'
  if (grepl("(pas|snowdepth)", var.name))
    leg.label <- 'cm'
  if (grepl("(tx90|tn10)", var.name))
    leg.label <- '%'
  if (grepl("(dd)", var.name))
    leg.label <- 'Degree days'
  if (grepl("(fd|cdd|cwd|su|gsl|id|trE|s30)", var.name))
    leg.label <- 'days'

  return(leg.label)
}

