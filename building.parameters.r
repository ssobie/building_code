##Contains functions to compute Building Code values from the BC
##Building Code and the American ASHRAE Handbook

library(ncdf4)
library(PCICt)
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

   R.dry <- 287.053 ##J K-1 kg-1
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
cdd.fxn  <- function(data,fac){tapply(data,fac, dd, tbase=10)}                            ##Cooling degree days
hdd.fxn  <- function(data,fac){tapply(-data,fac,dd, tbase=-10)}                           ##Heating degree days below 10
hdd.18.fxn  <- function(data,fac){tapply(-data,fac,dd, tbase=-18)}                           ##Heating degree days below 18

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
    return(list(ci=c(NA,NA,NA),
                rp=NA))
  } else {
    fx <- switch(var.name,
                 tasmax=max,
                 tasmin=min,
                 pr=max,
                 snow=max)

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

    ts.fit <- fevd(ts.to.fit,type='GEV')
    ts.old <- gev.fit(ts.to.fit,show=FALSE)

    ts.fit$results$par <- ts.old$mle
    names(ts.fit$results$par) <- c('location','scale','shape')
    ts.rp <- return.level(ts.fit,return.period=as.numeric(rperiod),make.plot=F)
    ts.rps <- c(NA,as.numeric(ts.rp),NA)

    tryCatch({
      ts.rps <- ci(ts.fit,alpha=0.05,type=c('return.level','parameter'),return.period=as.numeric(rperiod))
    }, warning=function(war) {
      message('Warning from confidence intervals')
    }, error=function(err) {
      print('Error in CI')
    }, finally={
    })

    ##This gives a 3 element vector with lower bound, rp value and upper bound
    rv <- as.numeric(ts.rps)
    rv <- list(ci=as.numeric(ts.rps),
               rp=as.numeric(ts.rp))

    if (var.name=='tasmin') {
      rv <- list(ci=-as.numeric(ts.rps)[c(3,2,1)],
                 rp=-as.numeric(ts.rp))
    }
    return(rv)
  }
}

##------------------------------------------------------------------
##Data subset

get.data.subset <- function(data,dates,interval) {

  bnds <- strsplit(interval,'-')[[1]]

  st <- head(grep(bnds[1],dates),1)
  en <- tail(grep(bnds[2],dates),1)

  data.subset <- data[st:en]
  dates.subset <- dates[st:en]

  rv <- list(data=data.subset,
             dates=dates.subset)
  return(rv)
}

get.variables.subset <- function(tas.data,tasmax.data,tasmin.data,
                                 pas.data,huss.data,pr.data,
                                 uas.data,vas.data,snd.data,dates,interval) {

  tas.mon <- get.data.subset(tas.data,dates,interval)
  tasmax.mon <- get.data.subset(tasmax.data,dates,interval) 
  tasmin.mon <- get.data.subset(tasmin.data,dates,interval)  

  pas.mon <- get.data.subset(pas.data,dates,interval) 
  huss.mon<- get.data.subset(huss.data,dates,interval) 
  pr.mon <- get.data.subset(pr.data,dates,interval) 

  uas.mon <- get.data.subset(uas.data,dates,interval) 
  vas.mon <- get.data.subset(vas.data,dates,interval) 

  snd.mon <- get.data.subset(vas.data,dates,interval) 

  rv <- list(tas=tas.mon,tasmax=tasmax.mon,tasmin=tasmin.mon,
             pas=pas.mon,huss=huss.mon,pr=pr.mon,
             uas=uas.mon,vas=vas.mon,snd=snd.mon)
  return(rv)  
}

ashrae.monthly.parameters <- function(data.subsets) {

  tas.sub <- data.subsets$tas                          
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin

  pas.sub <- data.subsets$pas 
  huss.sub <- data.subsets$huss  
  pr.sub <- data.subsets$pr 
                         
  dew.point.day <- dew.point.temp(pas.sub$data,huss.sub$data)
  wet.bulb.day  <- temp.wet.bulb(tas.sub$data,dew.point.day,pas.sub$data,huss.sub$data)
  wet.bulb.max  <- temp.wet.bulb(tasmax.sub$data,dew.point.day,pas.sub$data,huss.sub$data)
  wet.bulb.min  <- temp.wet.bulb(tasmin.sub$data,dew.point.day,pas.sub$data,huss.sub$data)

  avg.temp.mon <- mean.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%m')))
  sd.temp.mon <- sd.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%m')))

  hdd.sub.mon <- hdd.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%Y-%m')))
  hdd.sub <- tapply(hdd.sub.mon,as.factor(format(as.Date(paste(names(hdd.sub.mon),'-01',sep='')),'%m')),mean)
  cdd.sub.mon <- cdd.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%Y-%m')))
  cdd.sub <- tapply(cdd.sub.mon,as.factor(format(as.Date(paste(names(cdd.sub.mon),'-01',sep='')),'%m')),mean)
  cdd.24.sub.mon <- cdd.24.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%Y-%m')))
  cdd.24.sub <- tapply(cdd.24.sub.mon,as.factor(format(as.Date(paste(names(cdd.24.sub.mon),'-01',sep='')),'%m')),mean)

  avg.pr.mon <- mean.fxn(pr.sub$data,as.factor(format(pr.sub$dates,'%m')))
  max.pr.mon <- max.fxn(pr.sub$data,as.factor(format(pr.sub$dates,'%m')))
  min.pr.mon <- min.fxn(pr.sub$data,as.factor(format(pr.sub$dates,'%m')))
  sd.pr.mon <- sd.fxn(pr.sub$data,as.factor(format(pr.sub$dates,'%m')))
  
  wb.pctl.mon <- quant.fxn(wet.bulb.day,as.factor(format(tas.sub$dates,'%m')),pctl=0.996)
  wb.temp.mon <- conditional.mon.fxn(main=tas.sub$data,second=wet.bulb.day,fac=as.factor(format(tas.sub$dates,'%m')),pctl=0.996)

  db.pctl.mon <- quant.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%m')),pctl=0.996)
  db.temp.mon <- conditional.mon.fxn(main=wet.bulb.day,second=tas.sub$data,fac=as.factor(format(tas.sub$dates,'%m')),pctl=0.996)

  dtr.day <- tasmax.sub$data - tasmin.sub$data
  avg.dtr.mon <- mean.fxn(dtr.day,as.factor(format(tas.sub$dates,'%m')))
  dtr.temp.mon <- conditional.mon.fxn(main=dtr.day,second=tas.sub$data,fac=as.factor(format(tas.sub$dates,'%m')),pctl=0.95)

  wb.dtr.day <- wet.bulb.max - wet.bulb.min
  wb.dtr.mon <- conditional.mon.fxn(main=wb.dtr.day,second=tas.sub$data,fac=as.factor(format(tas.sub$dates,'%m')),pctl=0.95)
  print('Parameters')


  rv <- rbind(avg.temp.mon,sd.temp.mon,
              hdd.sub,cdd.sub,cdd.24.sub,
              avg.pr.mon,max.pr.mon,min.pr.mon,sd.pr.mon,
              wb.pctl.mon,wb.temp.mon,db.pctl.mon,db.temp.mon,
              avg.dtr.mon,dtr.temp.mon,wb.dtr.mon)

  return(rv)
}

ashrae.annual.wind.parameters <- function(data.subsets) {

  tas.sub <- data.subsets$tas                          
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin
  uas.sub <- data.subsets$uas 
  vas.sub <- data.subsets$vas 

  wspd <- sqrt(uas.sub$data^2 + vas.sub$data^2)
  wdir <- wind.dir(uas.sub$data,vas.sub$data)
  wspd.high.ann <- mean(quant.fxn(wspd,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  wspd.mcdb.low <- mean(conditional.mon.fxn(main=wspd,second=tasmin.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.004),na.rm=TRUE)
  db.mcws.low <- mean(conditional.mon.fxn(main=wdir,second=tasmin.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.004),na.rm=TRUE)

  wspd.mcdb.high <- mean(conditional.mon.fxn(main=wspd,second=tasmax.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  db.mcws.high <- mean(conditional.mon.fxn(main=wdir,second=tasmax.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)

  wspd.1.ann <- mean(quant.fxn(wspd,as.factor(format(tas.sub$dates,'%Y')),pctl=0.99),na.rm=TRUE)
  wspd.2.5.ann <- mean(quant.fxn(wspd,as.factor(format(tas.sub$dates,'%Y')),pctl=0.975),na.rm=TRUE)
  
  rv <- c(wspd.high.ann,wspd.mcdb.low,db.mcws.low)
  return(rv)

}

ashrae.annual.precip.parameters <- function(data.subsets) {

  tas.sub <- data.subsets$tas                          
  pr.sub <- data.subsets$pr                          

  pr.ann.avg <- mean(sum.fxn(pr.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.max <- max(sum.fxn(pr.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.min <- min(sum.fxn(pr.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.sd <- sd(sum.fxn(pr.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)

  rv <- c(pr.ann.avg,pr.ann.max,pr.ann.min,pr.ann.sd)  
  return(rv)
}

ashrae.annual.dry.temp.parameters <- function(data.subsets) {

  tas.sub <- data.subsets$tas                          
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin

  db.heat.ann <- mean(quant.fxn(tasmin.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.004),na.rm=TRUE)
  db.cool.ann <- mean(quant.fxn(tasmax.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)

  db.max.ann <- mean(max.fxn(tasmax.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  db.min.ann <- mean(min.fxn(tasmin.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)

  rv <- c(db.heat.ann,db.cool.ann,db.max.ann,db.min.ann)
  return(rv)
}

ashrae.return.periods <- function(data.subsets,rp) {

  tas.sub <- data.subsets$tas                          
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin
  pr.sub <- data.subsets$pr 

  tasmax.5.year <- calc.return.periods(tasmax.sub$data,tas.sub$dates,'tasmax',rp)
  tasmin.5.year <- calc.return.periods(tasmin.sub$data,tas.sub$dates,'tasmin',rp)
  pr.5.year <- calc.return.periods(pr.sub$data,tas.sub$dates,'pr',rp)  

  rv <- c(tasmax.5.year,tasmin.5.year,pr.5.year)
  return(rv)
}

ashrae.humidity.parameters <- function(data.subsets) {

  tas.sub <- data.subsets$tas                          
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin

  pas.sub <- data.subsets$pas 
  huss.sub <- data.subsets$huss  
  pr.sub <- data.subsets$pr 

  ##---
  dew.point.day <- dew.point.temp(pas.sub$data,huss.sub$data)
  wet.bulb.day  <- temp.wet.bulb(tas.sub$data,dew.point.day,pas.sub$data,huss.sub$data)

  dew.point.low <- mean(quant.fxn(dew.point.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.004),na.rm=TRUE)
  hum.ratio.low <- mean(quant.fxn(huss.sub$data*1000,as.factor(format(tas.sub$dates,'%Y')),pctl=0.004),na.rm=TRUE) ##In g/kg
  dew.mcdb.low <- mean(conditional.mon.fxn(main=tas.sub$data,second=dew.point.day,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.004),na.rm=TRUE)

  db.mcwb.ann <- mean(conditional.mon.fxn(main=wet.bulb.day,second=tasmax.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  wb.evap.ann <- mean(quant.fxn(wet.bulb.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  wb.mcdb.ann <- mean(conditional.mon.fxn(main=tas.sub$data,second=wet.bulb.day,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)

  dew.point.high <- mean(quant.fxn(dew.point.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  hum.ratio.high <- mean(quant.fxn(huss.sub$data*1000,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  dew.mcdb.high <- mean(conditional.mon.fxn(main=tas.sub$data,second=dew.point.day,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
 
  enth.day <- enthalpy(tas.sub$data,huss.sub$data)
  enth.ann <- mean(quant.fxn(enth.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  enth.mcdb.high <- mean(conditional.mon.fxn(main=tas.sub$data,second=enth.day,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE) 
  
  wb.max.ann <- mean(max.fxn(wet.bulb.day,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
      
  rv <- c(dew.point.low,hum.ratio.low,dew.mcdb.low,db.mcwb.ann,wb.evap.ann,wb.mcdb.ann,
          dew.point.high,hum.ratio.high,dew.mcdb.high,enth.ann,enth.mcdb.high,wb.max.ann)
  return(rv)

}

##------------------------------------------------------------------
##------------------------------------------------------------------
##BC Building Code

bc.return.periods <- function(data.subsets) {

  pr.sub <- data.subsets$pr 
  tas.sub <- data.subsets$tas
  snow.sub <- data.subsets$snd
  
  pr.rp.year <- calc.return.periods(pr.sub$data,tas.sub$dates,'pr',rp)  
  snow.rp.year <- calc.return.periods(snow.sub$data,tas.sub$dates,'snow',rp)  
  rv <- c(pr.rp.year,snow.rp.year)
  return(rv)

}

bc.quantiles <- function(data.subsets) {

  pr.sub <- data.subsets$pr 
  tas.sub <- data.subsets$tas
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin
  pas.sub <- data.subsets$pas
  huss.sub <- data.subsets$huss

  tasmin.025.mon <- quant.fxn(tasmin.sub,as.factor(format(tas.sub$dates,'%m')),pctl=0.025)
  tasmin.001.mon <- quant.fxn(tasmin.sub,as.factor(format(tas.sub$dates,'%m')),pctl=0.010)
  tasmax.975.mon <- quant.fxn(tasmin.sub,as.factor(format(tas.sub$dates,'%m')),pctl=0.975)

  dew.point.day <- dew.point.temp(pas.sub$data,huss.sub$data)
  wet.bulb.day  <- temp.wet.bulb(tas.sub$data,dew.point.day,pas.sub$data,huss.sub$data)
  wb.975.mon <- quant.fxn(tasmin.sub,as.factor(format(tas.sub$dates,'%m')),pctl=0.975)

  rv <- rbind(tasmin.025.mon,tasmin.001.mon,tasmax.975.mon,wb.975.mon)
  return(rv)
}

bc.mixed.parameters <- function(data.subsets) {

  pr.sub <- data.subsets$pr 
  tas.sub <- data.subsets$tas
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin

  hdd.18.ann <- mean(hdd.18.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.max <- mean(max.fxn(pr.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.avg <- mean(sum.fxn(pr.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)

  rv <- c(hdd.18.ann,pr.ann.max,pr.ann.avg)
  return(rv)
}

bc.wind.parameters <- function(data.subsets) {

  pr.sub <- data.subsets$pr 
  tas.sub <- data.subsets$tas
  uas.sub <- data.subsets$uas
  vas.sub <- data.subsets$vas
  huss.sub <- data.subsets$huss

  wind.rain.press.day <- wind.rain.press(uas,vas,pas,tas,huss)
  
  wind.press.5 <- calc.return.periods(wind.rain.press.day,tas.sub$dates,'pr',rp=5)    
  wind.press.10 <- calc.return.periods(wind.rain.press.day,tas.sub$dates,'pr',rp=10)    
  wind.press.50 <- calc.return.periods(wind.rain.press.day,tas.sub$dates,'pr',rp=50)    

  rv <- c(wind.press.5,wind.press.10,wind.press.50)
  return(rv)
}



##************************************************************************
##************************************************************************
##Load Data


read.dir <- '/storage/data/climate/downscale/BCCAQ2/CanRCM4/'

dates <- seq(from=as.Date('1950-01-01'),by='day',to=as.Date('2100-12-31'))
flag <- grep('*-02-29',dates)
dates <- dates[-flag]

load(paste(read.dir,'tas_NANAIMO_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
tas.data <- data - 273

load(paste(read.dir,'tasmax_NANAIMO_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
tasmax.data <- data

load(paste(read.dir,'tasmin_NANAIMO_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
tasmin.data <- data

load(paste(read.dir,'pr_NANAIMO_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
pr.data <- data
pr.data[pr.data < 0.1] <- 0

load(paste(read.dir,'ps_NANAIMO_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
pas.data <- data/1000

load(paste(read.dir,'huss_NANAIMO_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
huss.data <- data

load(paste(read.dir,'uas_NANAIMO_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
uas.data <- data

load(paste(read.dir,'vas_NANAIMO_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
vas.data <- data

load(paste(read.dir,'snd_NANAIMO_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
snd.data <- data


##-----------------------------------------------------------------------------
##Monthly Parameters
##Past

past.int <- '1971-2000'
past.var.subsets <- get.variables.subset(tas.data,tasmax.data,tasmin.data,
                                    pas.data,huss.data,pr.data,uas.data,vas.data,snd.data,
                                    dates,past.int)
past.mon.variables <- ashrae.monthly.parameters(past.var.subsets)
past.ann.variables <- ashrae.annual.parameters(past.var.subsets)
past.hum.variables <- ashrae.humidity.parameters(data.subsets)
past.dry.variables <- ashrae.annual.dry.temp.parameters(data.subsets)
past.wind.variables<- ashrae.annual.wind.parameters(data.subsets)
past.rp.variables  <- ashrae.return.periods(data.subsets,rp=5)


proj.int <- '2041-2070'
proj.var.subsets <- get.variables.subset(tas.data,tasmax.data,tasmin.data,
                                    pas.data,huss.data,pr.data,dates,proj.int)
proj.mon.variables <- ashrae.monthly.parameters(proj.var.subsets)

anom.mon.variables <- round(proj.mon.variables - past.mon.variables,1)
prct.mon.variables <- round((proj.mon.variables - past.mon.variables)/past.mon.variables * 100,1)

prct.mon.variables[is.infinite(prct.mon.variables) | is.nan(prct.mon.variables)] <- NA
prct.mon.variables[c(1,2),] <- NA

past.mon.variables <- round(past.mon.variables,1)
proj.mon.variables <- round(proj.mon.variables,1)

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
##BC Building Code



