##Contains functions to compute Building Code values from the BC
##Building Code and the American ASHRAE Handbook

library(ncdf4)
library(PCICt)

source('/storage/home/ssobie/code/repos/building_code/building.code.fxns.r')
source('/storage/home/ssobie/code/repos/building_code/convert.ashrae.to.metric.r')
source('/storage/home/ssobie/code/repos/building_code/gcm.bccaq.build.code.data.r',chdir=TRUE)

##------------------------------------------------------------------
##Nanaimo Variables
ashrae.vars <-  list(c(11.2,'degF'), ##Cold Dewpoint Temp 0.4% degF
                     c(9.8,'ratio'),  ##Cold Humidity Ratio 0.4% gr/lb
                     c(28.4,'degF'), ##MCDB Dewpoint 0.4% degF
                     c(63.9,'degF'), ##MCWB Dry Bulb 99.6% degF
                     c(62.6,'degF'), ##MCWB Dry Bulb 99.0% degF
                     c(61.4,'degF'), ##MCWB Dry Bulb 98.0% degF
                     c(65.2,'degF'), ##Evap Wet Bulb 99.6% degF
                     c(63.8,'degF'), ##Evap Wet Bulb 99.0% degF
                     c(76.8,'degF'), ##MCDB Wet Bulb 99.6% degF
                     c(73.9,'degF'), ##MCDB Wet Bulb 99.0% degF
                     c(60.8,'degF'), ##Warm Dewpoint 99.6% degF (Dehumidification)
                     c(59.6,'degF'), ##Warm Dewpoint 99.0% degF (Dehumidification)
                     c(79.9,'ratio'), ##Warm Humidity Ratio 99.6% gr/lb
                     c(76.3,'ratio'), ##Warm Humidity Ratio 99.0% gr/lb
                     c(68.9,'degF'),  ##MCDB Depoint 99.6% degF
                     c(67.4,'degF'),  ##MCDB Depoint 99.0% degF
                     c(30.1,'btu.lb'), ##Enthalpy 99.6% Btu/lb
                     c(77.0,'degF'),  ##MCDB Enthalpy 99.6% degF
                     c(75.2,'degF'), ##Annual Max Wet Bulb degF
                     c(23.4,'degF'), ##Heating Dry Bulb 0.4% degF
                     c(27.1,'degF'), ##Heating Dry Bulb 1.0% degF                     
                     c(80.3,'degF'), ##Cooling Dry Bulb 99.6% degF
                     c(76.5,'degF'), ##Cooling Dry Bulb 99.0% degF
                     c(87.5,'degF'), ##Annual Max Dry Bulb degF
                     c(20.1,'degF'), ##Annual Min Dry Bulb degF
                     c(46.4,'inch'), ##Annual Total Precip inches
                     c(62.4,'inch'), ##Annual Max Total Precip inches
                     c(28.1,'inch'), ##Annual Min Total Precip inches
                     c(8.0,'inch'),  ##Annual SD Total Precip inches
                     c(40.4,'mph'),  ##High Annual Wind Speed 99.6% mph
                     c(41.8,'degF'),  ##Windy MCDB Windy Speed 99.6% degF
                     c(6.8,'mph'),   ##Cold MCWS Dry Bulb 0.4% mph
                     c(300,'deg'),   ##Cold PCWD Dry Bulb 0.4% degrees from north
                     c(7.5,'mph'),   ##Warm MCWS Dry Bulb 99.6% mph
                     c(340,'deg'),   ##Warm PCWD Dry Bulb 99.6% degrees from north
                     c(21.9,'mph'),  ##Extreme Annual Wind Speed 95% mph
                     c(25.9,'mph'),  ##Extreme Annual Wind Speed 97.5% mph
                     c(30.2,'mph'),  ##Extreme Annual Wind Speed 99% mph
                     c(89.9,'degF'), ##5-Year Return Period Max Temp degF
                     c(17.0,'degF'), ##5-Year Return Period Min Temp degF
                     c(NA,'inch')) ##5-Year Return Period Daily Precip inches

ashrae.si.versions <- unlist(lapply(ashrae.vars,convert.to.si))

##Nanaimo Variables
bc.si.versions <-  c(-6.0, ##'degC'), ##Cold Month Design Temperature 2.5% degC
                     -8.0, ##'degC'),  ##Cold Month Design Temperature 1.0% degC
                     27,   ##'degC'), ##Warm Month Design Temperature 97.5% degC
                     19,   ##'degC'), ##Warm Month Wet Bulb Temperature 97.5% degC
                     3000, ##'DD'), ##Degree Days below 18 degC
                     NA,   ##'DD'), ##Degree Days below 18.3 degC
                     NA,   ##'DD'), ##Degree Days above 18.3 degC
                     NA,   ##'mm'), ##Annual Max Daily Precipitation mm
                     1050, ##'mm'), ##Annual Average Total Precipitation mm
                     91,   ##'mm'), ##50-Year return period Daily Precipitation
                     NA,  ##'Pa'),  ##5-Year return period Wind Load
                     390,  ##'Pa'))  ##10-Year return period Wind Load
                     500,  ##'Pa'),  ##50-Year return period Wind Load
                     NA,   ##'kPa'), ##20-Year Return Period Snow Load kPa
                     2.3)  ##,'kPa')) ##50-Year Return Period Snow Load kPa

##------------------------------------------------------------------

##Data subset

get.data.subset <- function(data,dates,interval) {

  bnds <- strsplit(interval,'-')[[1]]

  st <- head(grep(bnds[1],dates),1)
  en <- tail(grep(bnds[2],dates),1)
  if (length(en)==0) {
    en <- length(dates)
  }
  data.subset <- data[st:en]
  dates.subset <- dates[st:en]

  rv <- list(data=data.subset,
             dates=dates.subset)
  return(rv)
}

get.variables.subset <- function(tas.data,tasmax.data,tasmin.data,
                                 pas.data,huss.data,pr.data,
                                 uas.data,vas.data,snow.data,dates,interval) { ##snd.data,

  tas.mon <- get.data.subset(tas.data,dates,interval)
  tasmax.mon <- get.data.subset(tasmax.data,dates,interval) 
  tasmin.mon <- get.data.subset(tasmin.data,dates,interval)  

  pas.mon <- get.data.subset(pas.data,dates,interval) 
  huss.mon<- get.data.subset(huss.data,dates,interval) 
  pr.mon <- get.data.subset(pr.data,dates,interval) 

  uas.mon <- get.data.subset(uas.data,dates,interval) 
  vas.mon <- get.data.subset(vas.data,dates,interval) 

  snow.mon <- get.data.subset(snow.data,dates,interval) 

  dates.mon <- get.data.subset(dates,dates,interval)

  rv <- list(tas=tas.mon,tasmax=tasmax.mon,tasmin=tasmin.mon,
             pas=pas.mon,huss=huss.mon,pr=pr.mon,
             uas=uas.mon,vas=vas.mon,snow=snow.mon,dates=dates.mon) ##snd=snd.mon,
  return(rv)  
}

get.build.code.anomalies <- function(build.code.fxn,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets) {

  past.variables <- build.code.fxn(past.var.subsets)
  start.variables <- build.code.fxn(start.var.subsets)
  middle.variables <- build.code.fxn(middle.var.subsets)
  end.variables <- build.code.fxn(end.var.subsets)
  rx <- 2
  start.anom.variables <- round(start.variables - past.variables,rx)
  middle.anom.variables <- round(middle.variables - past.variables,rx)
  end.anom.variables <- round(end.variables - past.variables,rx)

  start.prct.variables <- round((start.variables - past.variables)/past.variables * 100,rx)
  middle.prct.variables <- round((middle.variables - past.variables)/past.variables * 100,rx)
  end.prct.variables <- round((end.variables - past.variables)/past.variables * 100,rx)

  start.prct.variables[is.infinite(start.prct.variables) | is.nan(start.prct.variables)] <- NA
  middle.prct.variables[is.infinite(middle.prct.variables) | is.nan(middle.prct.variables)] <- NA
  end.prct.variables[is.infinite(end.prct.variables) | is.nan(end.prct.variables)] <- NA

  rv <- list(past=past.variables,start.proj=start.variables,middle.proj=middle.variables,end.proj=end.variables,
             start.anom=start.anom.variables,middle.anom=middle.anom.variables,end.anom=end.anom.variables,
             start.prct=start.prct.variables,middle.prct=middle.prct.variables,end.prct=end.prct.variables)
  return(rv)
}



##------------------------------------------------------------------
##------------------------------------------------------------------
ashrae.rp.names <- rbind(c('5-Year Return Period Daily Max Temperature','degC'),   ##tasmax.5.year
                         c('5-Year Return Period Daily Min Temperature','degC'),   ##tasmin.5.year
                         c('5-Year Return Period Daily Total Precipitation','mm')) ##pr.5.year

ashrae.return.periods.5 <- function(data.subsets) {
  rp <- 5                        
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


ashrae.wind.names <- rbind(c('High Annual Wind Speed (99.6%)','m/s'),                             ##wspd.high
                           c('Mean Coincident Dry Bulb (Wind Speed 99.6%)','degC'),               ##mcdb.wspd.high
                           c('Cold Mean Coincident Wind Speed (Dry Bulb 0.4%)','m/s'),            ##wspd.high.ann
                           c('Cold Prevailing Coincident Wind Direction (Dry Bulb 0.4%)','deg'),  ##db.pcwd.low
                           c('Warm Mean Coincident Wind Speed (Dry Bulb 99.6%)','m/s'),           ##wspd.mcdb.ann
                           c('Warm Prevailing Coincident Wind Direction (Dry Bulb 99.6%)','deg'), ##db.pcwd.high
                           c('Extreme Annual Wind Speed (95%)','m/s'),                          ##wspd.5.ann
                           c('Extreme Annual Wind Speed (97.5%)','m/s'),                          ##wspd.2.5.ann
                           c('Extreme Annual Wind Speed (99%)','m/s'))                          ##wspd.1.ann

ashrae.annual.wind.parameters <- function(data.subsets) {

  tas.sub <- data.subsets$tas                          
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin
  uas.sub <- data.subsets$uas 
  vas.sub <- data.subsets$vas 

  wspd <- sqrt(uas.sub$data^2 + vas.sub$data^2)
  wdir <- wind.dir(uas.sub$data,vas.sub$data)
  wspd.high <- mean(quant.fxn(wspd,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  mcdb.wspd.high <- mean(conditional.mon.fxn(main=tas.sub$data,second=wspd,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)

  wspd.mcdb.low <- mean(conditional.mon.fxn(main=wspd,second=tasmin.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.004),na.rm=TRUE)
  db.pcwd.low <- mean(conditional.mon.fxn(main=wdir,second=tasmin.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.004),na.rm=TRUE)

  wspd.mcdb.high <- mean(conditional.mon.fxn(main=wspd,second=tasmax.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  db.pcwd.high <- mean(conditional.mon.fxn(main=wdir,second=tasmax.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)

  wspd.1.ann <- mean(quant.fxn(wspd,as.factor(format(tas.sub$dates,'%Y')),pctl=0.99),na.rm=TRUE)
  wspd.2.5.ann <- mean(quant.fxn(wspd,as.factor(format(tas.sub$dates,'%Y')),pctl=0.975),na.rm=TRUE)
  wspd.5.ann <- mean(quant.fxn(wspd,as.factor(format(tas.sub$dates,'%Y')),pctl=0.95),na.rm=TRUE)
  
  rv <- c(wspd.high,mcdb.wspd.high,wspd.mcdb.low,db.pcwd.low,wspd.mcdb.high,db.pcwd.high,wspd.5.ann,wspd.2.5.ann,wspd.1.ann)
  return(rv)
}

ashrae.precip.names <- rbind(c('Average Annual Total Precipitation','mm'),            ##pr.ann.avg
                             c('Maximum Annual Total Precipitation','mm'),            ##pr.ann.max
                             c('Minimum Annual Total Precipitation','mm'),            ##pr.ann.min
                             c('Standard Deviation Annual Total Precipitation','mm')) ##pr.ann.sd

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

ashrae.dry.temp.names <- rbind(c('Heating Dry Bulb Temperature 0.4%','degC'),   ##db.heat.ann
                               c('Heating Dry Bulb Temperature 1.0%','degC'),   ##db.heat.ann.1
                               c('Cooling Dry Bulb Temperature 99.6%','degC'),  ##db.cool.ann
                               c('Cooling Dry Bulb Temperature 99.0%','degC'),  ##db.cool.ann.1
                               c('Average Annual Max Dry Bulb Temperature','degC'),  ##db.max.ann
                               c('Average Annual Min Dry Bulb Temperature','degC'))  ##db.min.ann

ashrae.annual.dry.temp.parameters <- function(data.subsets) {

  tas.sub <- data.subsets$tas                          
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin

  db.heat.ann <- mean(quant.fxn(tasmin.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.004),na.rm=TRUE)
  db.heat.ann.1 <- mean(quant.fxn(tasmin.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.01),na.rm=TRUE)
  db.cool.ann <- mean(quant.fxn(tasmax.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  db.cool.ann.1 <- mean(quant.fxn(tasmax.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.99),na.rm=TRUE)

  db.max.ann <- mean(max.fxn(tasmax.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  db.min.ann <- mean(min.fxn(tasmin.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)

  rv <- c(db.heat.ann,db.heat.ann.1,db.cool.ann,db.cool.ann.1,db.max.ann,db.min.ann)
  return(rv)
}


ashrae.humidity.names <- rbind(c('Cold Dew Point Temperature 0.4%','degC'),            ##dew.point.low
                               c('Cold Humidity Ratio 0.4%','Ratio'),                       ##hum.ratio.low
                               c('Mean Coincident Dry Bulb (Dew Point 0.4%)','degC'),  ##dew.mcdb.low
                               c('Mean Coincident Wet Bulb (Dry Bulb 99.6%)','degC'),  ##db.mcwb.ann.996
                               c('Mean Coincident Wet Bulb (Dry Bulb 99.0%)','degC'),  ##db.mcwb.ann.99
                               c('Mean Coincident Wet Bulb (Dry Bulb 98.0%)','degC'),  ##db.mcwb.ann.98
                               c('Evaporation Wet Bulb 99.6%','degC'),                 ##wb.evap.ann.996
                               c('Evaporation Wet Bulb 99.0%','degC'),                 ##wb.evap.ann.99
                               c('Mean Coincident Dry Bulb (Wet Bulb 99.6%)','degC'),  ##wb.mcdb.ann.996
                               c('Mean Coincident Dry Bulb (Wet Bulb 99.0%)','degC'),  ##wb.mcdb.ann.99
                               c('Warm Dew Point Temperature 99.6%','degC'),           ##dew.point.996
                               c('Warm Dew Point Temperature 99.0%','degC'),           ##dew.point.990
                               c('Warm Humidity Ratio 99.6%','Ratio'),                 ##hum.ratio.996
                               c('Warm Humidity Ratio 99.0%','Ratio'),                 ##hum.ratio.990
                               c('Mean Coincident Dry Bulb (Dew Point 99.6%)','degC'), ##dew.mcdb.996
                               c('Mean Coincident Dry Bulb (Dew Point 99.0%)','degC'), ##dew.mcdb.990
                               c('Enthalpy 99.6%','kJ/kg'),                           ##enth.ann
                               c('Mean Coincident Dry Bulb (Enthalpy 99.6%)','degC'),  ##enth.mcdb.high
                               c('Annual Maximum Wet Bulb','degC'))                    ##wb.max.ann       

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

  db.mcwb.ann.996 <- mean(conditional.mon.fxn(main=wet.bulb.day,second=tasmax.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  db.mcwb.ann.99 <- mean(conditional.mon.fxn(main=wet.bulb.day,second=tasmax.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.99),na.rm=TRUE)
  db.mcwb.ann.98 <- mean(conditional.mon.fxn(main=wet.bulb.day,second=tasmax.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.98),na.rm=TRUE)
  wb.evap.ann.996 <- mean(quant.fxn(wet.bulb.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  wb.evap.ann.99 <- mean(quant.fxn(wet.bulb.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.99),na.rm=TRUE)
  wb.mcdb.ann.996 <- mean(conditional.mon.fxn(main=tas.sub$data,second=wet.bulb.day,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  wb.mcdb.ann.99 <- mean(conditional.mon.fxn(main=tas.sub$data,second=wet.bulb.day,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.99),na.rm=TRUE)

  dew.point.996 <- mean(quant.fxn(dew.point.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  dew.point.990 <- mean(quant.fxn(dew.point.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.99),na.rm=TRUE)
  hum.ratio.996 <- mean(quant.fxn(huss.sub$data*1000,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  hum.ratio.990 <- mean(quant.fxn(huss.sub$data*1000,as.factor(format(tas.sub$dates,'%Y')),pctl=0.990),na.rm=TRUE)
  dew.mcdb.996 <- mean(conditional.mon.fxn(main=tas.sub$data,second=dew.point.day,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  dew.mcdb.990 <- mean(conditional.mon.fxn(main=tas.sub$data,second=dew.point.day,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.990),na.rm=TRUE)
 
  enth.day <- enthalpy(tas.sub$data,huss.sub$data)
  enth.ann <- mean(quant.fxn(enth.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)/1000 ##to convert to kJ
  enth.mcdb.high <- mean(conditional.mon.fxn(main=tas.sub$data,second=enth.day,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE) 
  
  wb.max.ann <- mean(max.fxn(wet.bulb.day,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  rv <- c(dew.point.low,hum.ratio.low,dew.mcdb.low,db.mcwb.ann.996,db.mcwb.ann.99,db.mcwb.ann.98,wb.evap.ann.996,wb.evap.ann.99,wb.mcdb.ann.996,wb.mcdb.ann.99,
          dew.point.996,dew.point.990,hum.ratio.996,hum.ratio.990,dew.mcdb.996,dew.mcdb.990,enth.ann,enth.mcdb.high,wb.max.ann)
  return(rv)

}

##------------------------------------------------------------------
##------------------------------------------------------------------
##BC Building Code

bc.quantile.names <- rbind(c('Cool Month Design Temperature 2.5%','degC'), ##tasmin.025.mon
                           c('Cold Month Design Temperature 1.0%','degC'), ##tasmin.001.mon
                           c('Warm Month Design Temperature 97.5%','degC'), ##tasmax.975.mon
                           c('Warm Month Design Wet Bulb Temp 97.5%','degC')) ##wb.975.mon

bc.quantile.parameters <- function(data.subsets) {

  pr.sub <- data.subsets$pr 
  tas.sub <- data.subsets$tas
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin
  pas.sub <- data.subsets$pas
  huss.sub <- data.subsets$huss

  tasmin.025.mon <- mean(quant.fxn(tasmin.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.025),na.rm=TRUE)
  tasmin.001.mon <- mean(quant.fxn(tasmin.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.010),na.rm=TRUE)
  tasmax.975.mon <- mean(quant.fxn(tasmax.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.975),na.rm=TRUE)

  dew.point.day <- dew.point.temp(pas.sub$data,huss.sub$data)
  wet.bulb.day  <- temp.wet.bulb(tas.sub$data,dew.point.day,pas.sub$data,huss.sub$data)
  wb.975.mon <- mean(quant.fxn(wet.bulb.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.975),na.rm=TRUE)

  rv <- rbind(tasmin.025.mon,tasmin.001.mon,tasmax.975.mon,wb.975.mon)
  return(rv)
}

bc.mixed.names <- rbind(c('Heating Degree Days (18 degC)','DD'), ##hdd.18.ann
                        c('Heating Degree Days (18.3 degC)','DD'), ##hdd.18.3.ann
                        c('Cooling Degree Days (18.3 degC)','DD'), ##cdd.18.3.ann
                        c('Annual Maximum Daily Precipitation','mm'), ##pr.ann.max
                        c('Annual Average Total Precipitation','mm'),
                        c('50-Year Return Period Daily Total Precipitation','mm')) ##pr.rp.year


bc.mixed.parameters <- function(data.subsets) {

  pr.sub <- data.subsets$pr 
  tas.sub <- data.subsets$tas
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin

  hdd.18.ann <- mean(hdd.18.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  hdd.18.3.ann <- mean(hdd.18.3.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  cdd.18.3.ann <- mean(cdd.18.3.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.max <- mean(max.fxn(pr.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.avg <- mean(sum.fxn(pr.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  pr.rp.year <- calc.return.periods(pr.sub$data,tas.sub$dates,'pr',rp=50)  

  rv <- c(hdd.18.ann,hdd.18.3.ann,cdd.18.3.ann,pr.ann.max,pr.ann.avg,pr.rp.year)
  return(rv)
}

bc.snow.names <- rbind(c('20-Year Return Period Daily Snow Load','kPa'), ##pr.rp.year
                       c('50-Year Return Period Daily Snow Load','kPa')) ##pr.rp.year
bc.snow.parameters <- function(data.subsets) {
  snow.sub <- data.subsets$snow
  snow.rp20.year <- calc.return.periods(snow.sub$data*100,snow.sub$dates,'snd',rp=20)*3.5/100
  snow.rp50.year <- calc.return.periods(snow.sub$data*100,snow.sub$dates,'snd',rp=50)*3.5/100
  ##snow.rp.year <- snow.rp.year*3.5 ##To convert to snow load assuming snow density of 3.5kN/m^3  
  return(c(snow.rp20.year,snow.rp50.year))
}


bc.wind.names <- rbind(c('5-Year Wind Pressure Return Period','Pa'), ##wind.press.5
                       c('10-Year Wind Pressure Return Period','Pa'), ##wind.press.10
                       c('50-Year Wind Pressure Return Period','Pa')) ##wind.press.50

bc.wind.parameters <- function(data.subsets) {

  pr.sub <- data.subsets$pr 
  tas.sub <- data.subsets$tas
  uas.sub <- data.subsets$uas
  vas.sub <- data.subsets$vas
  huss.sub <- data.subsets$huss
  pas.sub <- data.subsets$pas

  wind.rain.press.day <- wind.rain.press(uas.sub$data,vas.sub$data,pas.sub$data,tas.sub$data,huss.sub$data)
  wind.press.5 <- calc.return.periods(wind.rain.press.day,tas.sub$dates,'pr',rp=5)    
  wind.press.10 <- calc.return.periods(wind.rain.press.day,tas.sub$dates,'pr',rp=10)    
  wind.press.50 <- calc.return.periods(wind.rain.press.day,tas.sub$dates,'pr',rp=50)    

  rv <- c(wind.press.5,wind.press.10,wind.press.50)*1000
  return(rv)
}

##------------------------------------------------------------------
##------------------------------------------------------------------

get.var.title <- function(var.name) {

  var.title <- switch(var.name,
                      mon.cool='TASMIN 2.5% degC',
                      mon.cold='TASMIN 1.0% degC',
                      dry.warm='TASMAX 97.5% degC',
                      wet.warm='WET BULB 97.5% degC',
                      avg.temp.mon='Average Temp degC',
                      sd.temp.mon='St.Dev. Temp degC',
                      hdd.sub='Heating Degree Days (10 degC)',    
                      cdd.sub='Cooling Degree Days (10 degC)',   
                      cdd.24.sub='Cooling Degree Days (24 degC)',    
                      avg.pr.mon='Avg Precipitation', 
                      max.pr.mon='Max Precipitation',
                      min.pr.mon='Min Precipitation',
                      sd.pr.mon='St.Dev. Precipitation',  
                      wb.pctl.mon='Wet Bulb 99.6% degC',
                      wb.temp.mon='Mean Conditional Dry Bulb 99.6% degC', 
                      db.pctl.mon='Dry Bulb 99.6% degC',
                      db.temp.mon='Mean Conditional Wet Bulb 99.6% degC', 
                      avg.dtr.mon='Diurnal Temperature Range degC', 
                      dtr.temp.mon='Mean Conditional Dry Bulb Range 95% degC', 
                      wb.dtr.mon= 'Mean Conditional Wet Bulb Range 95% degC',
                      tasmin.025.mon='Cold Design Temperature 2.5% degC',       
                      tasmin.001.mon='Cold Design Temperature 1.0% degC',
                      tasmax.975.mon='Warm Design Temperature 97.5% degC',
                      wb.975.mon='Warm Wet Bulb Design Temperature 97.5% degC',
                      dew.point.low='Cold Dew Point Temperature 0.4% degC',
                      hum.ratio.low='Cold Humidity Ratio 0.4%',
                      dew.mcdb.low='Mean Conditional Dry Bulb (Dew Point 0.4%) degC',
                      db.mcwb.ann='Mean Conditional Wet Bulb (Dry Bulb 0.4%) degC',
                      wb.evap.ann='Evaporation Wet Bulb ')

  return(var.title)              
}


monthly.table <- function(var.title,var.units,data) {
  data.mon <- cbind(month.abb,data)              
  data.bottom <- c('Month',paste('Past (',var.units,')',sep=''),
                 paste('2020s Projection (',var.units,')',sep=''),
                 paste('2020s Change (',var.units,')',sep=''),
                 paste('2020s Percent Change (%)',sep=''),                       
                 paste('2050s Projection (',var.units,')',sep=''),               
                 paste('2050s Change (',var.units,')',sep=''),
                 paste('2050s Percent Change (%)',sep=''), 
                 paste('2080s Projection (',var.units,')',sep=''),               
                 paste('2080s Change (',var.units,')',sep=''),
                 paste('2080s Percent Change (%)',sep=''))
  data.top <- c(var.title,rep(' ',length(data.bottom)-1))                 
  data.header <- rbind(data.top,data.bottom)
  var.result <- rbind(data.header,data.mon)
  return(var.result)
}

annual.table <- function(var.title,var.units,data) {

  data.header  <- c(var.title,
                 paste('Handbook Past (',var.units,')',sep=''),
                 paste('Model Past (',var.units,')',sep=''),
                 paste('2020s Projection (',var.units,')',sep=''),
                 paste('2020s Change (',var.units,')',sep=''),
                 paste('2020s Percent Change (%)',sep=''),                       
                 paste('2050s Projection (',var.units,')',sep=''),               
                 paste('2050s Change (',var.units,')',sep=''),
                 paste('2050s Percent Change (%)',sep=''), 
                 paste('2080s Projection (',var.units,')',sep=''),               
                 paste('2080s Change (',var.units,')',sep=''),
                 paste('2080s Percent Change (%)',sep=''))
  var.result <- rbind(data.header,c('',round(data,2)))

  if (var.units=='degC') {
     var.result[2,c(6,9,12)] <- NA
  }
  return(var.result)
}

write.annual.table <- function(data.list,data.names,write.file) {

   data.matrix <- cbind(data.list$site,data.list$past,
                        data.list$start.proj,data.list$start.anom,data.list$start.prct,
                        data.list$middle.proj,data.list$middle.anom,data.list$middle.prct,
                        data.list$end.proj,data.list$end.anom,data.list$end.prct)

   data.len <- nrow(data.matrix)
   data.vars <- c()
   for (i in 1:data.len) {
       data.formatted <- annual.table(var.title=data.names[i,1],
                                        var.units=data.names[i,2],
                                        data=data.matrix[i,])
       data.vars <- rbind(data.vars,data.formatted)
   }

   my.writedir <- '/storage/data/projects/rci/building_code/'
   write.table(data.vars,file=write.file,
               sep=',',quote=FALSE,col.name=FALSE,row.name=FALSE)
}

write.combined.table <- function(means.list,tenth.list,ninetieth.list,data.names,write.file) {

   range.matrix <- cbind(means.list$site,
                         means.list$past,paste(tenth.list$past,ninetieth.list$past,sep=' to '),
                         means.list$start.proj,paste(tenth.list$start.proj,ninetieth.list$start.proj,sep=' to '),
                         means.list$start.anom,paste(tenth.list$start.anom,ninetieth.list$start.anom,sep=' to '),
                         means.list$start.prct,paste(tenth.list$start.prct,ninetieth.list$start.prct,sep=' to '),
                         means.list$middle.proj,paste(tenth.list$middle.proj,ninetieth.list$middle.proj,sep=' to '),
                         means.list$middle.anom,paste(tenth.list$middle.anom,ninetieth.list$middle.anom,sep=' to '),
                         means.list$middle.prct,paste(tenth.list$middle.prct,ninetieth.list$middle.prct,sep=' to '),
                         means.list$end.proj,paste(tenth.list$end.proj,ninetieth.list$end.proj,sep=' to '),
                         means.list$end.anom,paste(tenth.list$end.anom,ninetieth.list$end.anom,sep=' to '),
                         means.list$end.prct,paste(tenth.list$end.prct,ninetieth.list$end.prct,sep=' to '))

    data.header  <- rbind(c('Handbook Past','Model Past',' ','2020s Projection',' ','2020s Change',' ','2020s Percent Change',' ',
                                                                           '2050s Projection',' ','2050s Change',' ','2050s Percent Change',' ',
                                                                           '2080s Projection',' ','2080s Change',' ','2080s Percent Change',' '),
                        c('Site','Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%',
                                                             'Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%',
                                                             'Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%'))
   title.col <- cbind(c('Building Code','Variable',data.names[,1]),c('Variable','Units',data.names[,2]))

   data.table <- rbind(data.header,range.matrix)
   data.out <- cbind(title.col,data.table)                 
   
   ts.ix <- grep('degC',data.out[,2])
   pc.ix <- grep('Percent',data.out[1,])
   data.out[ts.ix,c(pc.ix,pc.ix+1)] <- NA

   my.writedir <- '/storage/data/projects/rci/building_code/'
   write.table(data.out,file=write.file,
               sep=',',quote=FALSE,col.name=FALSE,row.name=FALSE)
}


##************************************************************************
##************************************************************************

model.list <- c('ACCESS1-0',
              'CanESM2',
              'CNRM-CM5',
              'CSIRO-Mk3-6-0',
              'GFDL-ESM2G',
              'HadGEM2-CC',
              'HadGEM2-ES',
              'inmcm4',
              'MIROC5',
              'MRI-CGCM3')

rcp26.list <- c('CanESM2',
                'CNRM-CM5',
                'CSIRO-Mk3-6-0',
                'GFDL-ESM2G',
                'HadGEM2-ES',
                'MIROC5',
                'MRI-CGCM3')

scen.list <- 'rcp85' ##c('rcp26','rcp45','rcp85')

data.dir <- '/storage/data/climate/downscale/CMIP5/building_code/data_files/'

my.writedir <- '/storage/data/projects/rci/building_code/'

##Nanaimo Hospital Coordinates
lon.c <- -123.969357
lat.c <- 49.184737

for (scenario in scen.list) {
    gcm.list <- model.list
    if (scenario=='rcp26') {
       gcm.list <- rcp26.list
    }
    ashrae.matrix <- vector(mode='list',length=length(gcm.list))
    bc.matrix <- vector(mode='list',length=length(gcm.list))
    
    for (g in seq_along(gcm.list)) {
        gcm <- gcm.list[g]
        print(gcm)
        ##BCCAQ Data
        bccaq.file <- paste0(data.dir,'nanaimo_hospital_bccaq_',gcm,'_',scenario,'_pr_tx_tn.RData')
        if (!file.exists(bccaq.file)) {
           bccaq.data <- gather.bccaq.data(gcm,lon.c,lat.c,scenario)        
           save(bccaq.data,file=bccaq.file)
        } else {
           load(bccaq.file)
        }
        if (scenario=='rcp85') {
          snow.file <-  paste0(data.dir,'nanaimo_hospital_bccaq_',gcm,'_',scenario,'_snow.RData')
          if (!file.exists(snow.file)) {
             snow.data <- gather.snow.data(gcm,lon.c,lat.c,scenario)        
             save(snow.data,file=snow.file)
          } else {
             load(snow.file)
          }
        }
        tas.data <- bccaq.data$tas
        tasmax.data <- bccaq.data$tasmax
        tasmin.data <- bccaq.data$tasmin
        pr.data <- bccaq.data$pr        
        snow.data <- snow.data$snow

        ##GCM other variables
        gcm.file <- paste0(data.dir,'nanaimo_hospital_gcm_',gcm,'_',scenario,'_huss_psl_uas_vas.RData')
        if (!file.exists(gcm.file)) {
           gcm.data <- gather.gcm.data(gcm,lon.c,lat.c,scenario)        
           save(gcm.data,file=gcm.file)
        } else {
           load(gcm.file)
        }

        huss.data <- gcm.data$huss
        psl.data  <- gcm.data$psl
        uas.data  <- gcm.data$uas
        vas.data  <- gcm.data$vas       

        ##-----------------------------------------------------------------------------
        ##Monthly Parameters
        ##Past
        past.int <- '1971-2000'
        past.var.subsets <- get.variables.subset(tas.data$data,tasmax.data$data,tasmin.data$data,
                                    psl.data$data,huss.data$data,pr.data$data,uas.data$data,vas.data$data,snow.data$data,
                                    tasmax.data$time,past.int)

        start.int <- '2011-2040'
        start.var.subsets <- get.variables.subset(tas.data$data,tasmax.data$data,tasmin.data$data,
                                    psl.data$data,huss.data$data,pr.data$data,uas.data$data,vas.data$data,snow.data$data,
                                    tasmax.data$time,start.int)
                                    
        middle.int <- '2041-2070'
        middle.var.subsets <- get.variables.subset(tas.data$data,tasmax.data$data,tasmin.data$data,
                                    psl.data$data,huss.data$data,pr.data$data,uas.data$data,vas.data$data,snow.data$data,
                                    tasmax.data$time,middle.int)

        end.int <- '2071-2100'
        end.var.subsets <- get.variables.subset(tas.data$data,tasmax.data$data,tasmin.data$data,
                                    psl.data$data,huss.data$data,pr.data$data,uas.data$data,vas.data$data,snow.data$data,
                                    tasmax.data$time,end.int)

        ##------------------------------------------------------------------------------------------------
        ##------------------------------------------------------------------------------------------------
        ##ASHRAE Parameters
if (1==0) {
        ashrae.humidity <- get.build.code.anomalies(ashrae.humidity.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)

        ashrae.dry.temp <- get.build.code.anomalies(ashrae.annual.dry.temp.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)

        ashrae.precip   <- get.build.code.anomalies(ashrae.annual.precip.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)

        ashrae.wind     <- get.build.code.anomalies(ashrae.annual.wind.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)

        ashrae.rp       <- get.build.code.anomalies(ashrae.return.periods.5,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)

        ashrae.all.parameters <- mapply(c,ashrae.humidity,ashrae.dry.temp,ashrae.precip,ashrae.wind,ashrae.rp,SIMPLIFY=FALSE)
        ashrae.all.names <- rbind(ashrae.humidity.names,
                          ashrae.dry.temp.names,
                          ashrae.precip.names,
                          ashrae.wind.names,
                          ashrae.rp.names)

        ashrae.all.parameters <- c(list(site=round(ashrae.si.versions,1)),ashrae.all.parameters)
        if (g==1) {
           ashrae.matrix <- ashrae.all.parameters
        } else {
           ashrae.matrix <- mapply(rbind,ashrae.matrix,ashrae.all.parameters,SIMPLIFY=FALSE)
        }
}        
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
        ##BC Parameters


        bc.quantiles <- get.build.code.anomalies(bc.quantile.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
        bc.mixed <- get.build.code.anomalies(bc.mixed.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
        print('SNOW')
        bc.snow <- get.build.code.anomalies(bc.snow.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
        if (gcm=='HadGEM2-ES') {
          bc.snow <- lapply(bc.snow,function(x){x*NA})
        }
        bc.wind <- get.build.code.anomalies(bc.wind.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
        
        bc.all.parameters <- mapply(c,bc.quantiles,bc.mixed,bc.wind,bc.snow,SIMPLIFY=FALSE)
        bc.all.names <- rbind(bc.quantile.names,
                          bc.mixed.names,
                          bc.wind.names,
                          bc.snow.names)

        bc.all.parameters <- c(list(site=round(bc.si.versions,2)),bc.all.parameters)

        if (g==1) {
           bc.matrix <- bc.all.parameters
        } else {
           bc.matrix <- mapply(rbind,bc.matrix,bc.all.parameters,SIMPLIFY=FALSE)        
        }
     }
if (1==0) {
        ashrae.means <- lapply(lapply(ashrae.matrix,function(x){apply(x,2,mean)}),round,1)
        ashrae.means$start.prct <- round(ashrae.means$start.anom/ashrae.means$past*100,1)
        ashrae.means$middle.prct <- round(ashrae.means$middle.anom/ashrae.means$past*100,1)
        ashrae.means$end.prct <- round(ashrae.means$end.anom/ashrae.means$past*100,1)

##        write.file <- paste(my.writedir,'ashrae.gcm.bccaq.',scenario,'.ens.means.all.variables.nanaimo.building.code.csv',sep='')
##        write.annual.table(ashrae.means,ashrae.all.names,write.file)

        ashrae.10th <- lapply(lapply(ashrae.matrix,function(x){apply(x,2,quantile,0.1,na.rm=T)}),round,1)
        ashrae.10th$start.prct <- round(ashrae.10th$start.anom/ashrae.10th$past*100,1)
        ashrae.10th$middle.prct <- round(ashrae.10th$middle.anom/ashrae.10th$past*100,1)
        ashrae.10th$end.prct <- round(ashrae.10th$end.anom/ashrae.10th$past*100,1)

##        write.10th <- paste(my.writedir,'ashrae.gcm.bccaq.',scenario,'.ens.10th.percentile.all.variables.nanaimo.building.code.csv',sep='')
##        write.annual.table(ashrae.10th,ashrae.all.names,write.10th)

        ashrae.90th <- lapply(lapply(ashrae.matrix,function(x){apply(x,2,quantile,0.9,na.rm=T)}),round,1)
        ashrae.90th$start.prct <- round(ashrae.90th$start.anom/ashrae.90th$past*100,1)
        ashrae.90th$middle.prct <- round(ashrae.90th$middle.anom/ashrae.90th$past*100,1)
        ashrae.90th$end.prct <- round(ashrae.90th$end.anom/ashrae.90th$past*100,1)

##        write.90th <- paste(my.writedir,'ashrae.gcm.bccaq.',scenario,'.ens.90th.percentile.all.variables.nanaimo.building.code.csv',sep='')
##        write.annual.table(ashrae.90th,ashrae.all.names,write.90th)

        write.combined.file <- paste(my.writedir,'ashrae.gcm.bccaq.',scenario,'.ens.combined.all.variables.nanaimo.building.code.csv',sep='')
        write.combined.table(ashrae.means,ashrae.10th,ashrae.90th,ashrae.all.names,write.combined.file)
}
        bc.means <- lapply(lapply(bc.matrix,function(x){apply(x,2,mean,na.rm=T)}),round,1)
        bc.means$start.prct <- round(bc.means$start.anom/bc.means$past*100,1)
        bc.means$middle.prct <- round(bc.means$middle.anom/bc.means$past*100,1)
        bc.means$end.prct <- round(bc.means$end.anom/bc.means$past*100,1)
##        write.file <- paste(my.writedir,'bc.gcm.bccaq.',scenario,'.ens.means.all.variables.nanaimo.building.code.csv',sep='')
##        write.annual.table(bc.means,bc.all.names,write.file)

        bc.10th <- lapply(lapply(bc.matrix,function(x){apply(x,2,quantile,0.1,na.rm=T)}),round,1)
        bc.10th$start.prct <- round(bc.10th$start.anom/bc.10th$past*100,1)
        bc.10th$middle.prct <- round(bc.10th$middle.anom/bc.10th$past*100,1)
        bc.10th$end.prct <- round(bc.10th$end.anom/bc.10th$past*100,1)

##        write.10th <- paste(my.writedir,'bc.gcm.bccaq.',scenario,'.ens.10th.percentile.all.variables.nanaimo.building.code.csv',sep='')
##        write.annual.table(bc.10th,bc.all.names,write.10th)

        bc.90th <- lapply(lapply(bc.matrix,function(x){apply(x,2,quantile,0.9,na.rm=T)}),round,1)
        bc.90th$start.prct <- round(bc.90th$start.anom/bc.90th$past*100,1)
        bc.90th$middle.prct <- round(bc.90th$middle.anom/bc.90th$past*100,1)
        bc.90th$end.prct <- round(bc.90th$end.anom/bc.90th$past*100,1)

##        write.90th <- paste(my.writedir,'bc.gcm.bccaq.',scenario,'.ens.90th.percentile.all.variables.nanaimo.building.code.csv',sep='')
##        write.annual.table(bc.90th,bc.all.names,write.90th)

        write.combined.file <- paste(my.writedir,'bc.gcm.bccaq.',scenario,'.ens.combined.all.variables.nanaimo.building.code.csv',sep='')
        write.combined.table(bc.means,bc.10th,bc.90th,bc.all.names,write.combined.file)


}



