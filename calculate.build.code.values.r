##Contains functions to compute Building Code values from the BC
##Building Code and the American ASHRAE Handbook

library(ncdf4)
library(PCICt)

source('/storage/home/ssobie/code/repos/building_code/building.code.fxns.r')
source('/storage/home/ssobie/code/repos/building_code/convert.ashrae.to.metric.r')
source('/storage/home/ssobie/code/repos/building_code/gcm.bccaq.build.code.data.r',chdir=TRUE)
source('/storage/home/ssobie/code/repos/building_code/formatted.building.code.table.r',chdir=TRUE)

##------------------------------------------------------------------
##Nanaimo ASHRAE Variables and Duncan BCBC Variables 
##for Cowichan Hospital
##source('/storage/home/ssobie/code/repos/building_code/cowichan.hospital.build.observations.r')
source('/storage/home/ssobie/code/repos/building_code/granville.41st.build.observations.r')

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
             dates=as.Date(dates.subset))
  return(rv)
}

get.variables.subset <- function(tas.data,tasmax.data,tasmin.data,
                                 pas.data,huss.data,pr.data,
                                 uas.data,vas.data,snow.data,interval) { ##snd.data,

  tas.mon <- get.data.subset(tas.data$data,tas.data$time,interval)
  tasmax.mon <- get.data.subset(tasmax.data$data,tasmax.data$time,interval) 
  tasmin.mon <- get.data.subset(tasmin.data$data,tasmin.data$time,interval)  

  pas.mon <- get.data.subset(pas.data$data,pas.data$time,interval) 
  huss.mon<- get.data.subset(huss.data$data,huss.data$time,interval) 
  pr.mon <- get.data.subset(pr.data$data,pr.data$time,interval) 

  uas.mon <- get.data.subset(uas.data$data,uas.data$time,interval) 
  vas.mon <- get.data.subset(vas.data$data,vas.data$time,interval) 

  snow.mon <- get.data.subset(snow.data$data,snow.data$time,interval) 

  ##dates.mon <- get.data.subset(dates,dates,interval)

  rv <- list(tas=tas.mon,tasmax=tasmax.mon,tasmin=tasmin.mon,
             pas=pas.mon,huss=huss.mon,pr=pr.mon,
             uas=uas.mon,vas=vas.mon,snow=snow.mon) ##snd=snd.mon,
  return(rv)  
}

fix.hadley <- function(data) {
  new.data <- data$data[1:10440]
  new.dates <- data$dates[1:10440]
  rv <- list(data=new.data,dates=new.dates)
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
  pr.5.year <- calc.return.periods(pr.sub$data,pr.sub$dates,'pr',rp)  

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
  wspd.high <- mean(quant.fxn(wspd,as.factor(format(uas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  mcdb.wspd.high <- mean(conditional.mon.fxn(main=tas.sub$data,second=wspd,fac=as.factor(format(uas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)

  wspd.mcdb.low <- mean(conditional.mon.fxn(main=wspd,second=tasmin.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.004),na.rm=TRUE)
  db.pcwd.low <- mean(conditional.mon.fxn(main=wdir,second=tasmin.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.004),na.rm=TRUE)

  wspd.mcdb.high <- mean(conditional.mon.fxn(main=wspd,second=tasmax.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  db.pcwd.high <- mean(conditional.mon.fxn(main=wdir,second=tasmax.sub$data,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)

  wspd.1.ann <- mean(quant.fxn(wspd,as.factor(format(uas.sub$dates,'%Y')),pctl=0.99),na.rm=TRUE)
  wspd.2.5.ann <- mean(quant.fxn(wspd,as.factor(format(uas.sub$dates,'%Y')),pctl=0.975),na.rm=TRUE)
  wspd.5.ann <- mean(quant.fxn(wspd,as.factor(format(uas.sub$dates,'%Y')),pctl=0.95),na.rm=TRUE)
  
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

  pr.ann.avg <- mean(sum.fxn(pr.sub$data,as.factor(format(pr.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.max <- max(sum.fxn(pr.sub$data,as.factor(format(pr.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.min <- min(sum.fxn(pr.sub$data,as.factor(format(pr.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.sd <- sd(sum.fxn(pr.sub$data,as.factor(format(pr.sub$dates,'%Y'))),na.rm=TRUE)

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

bc.quantile.names <- rbind(c('tasmin.050.mon','Cool Month Design Temperature 5.0%','degC'), ##tasmin.050.mon
                           c('tasmin.025.mon','Cool Month Design Temperature 2.5%','degC'), ##tasmin.025.mon
                           c('tasmin.010.mon','Cold Month Design Temperature 1.0%','degC'), ##tasmin.010.mon
                           c('tasmin.004.mon','Cold Month Design Temperature 0.4%','degC'), ##tasmin.004.mon
                           c('tasmax.996.mon','Warm Month Design Temperature 99.6%','degC'), ##tasmax.996.mon
                           c('tasmax.990.mon','Warm Month Design Temperature 99.0%','degC'), ##tasmax.990.mon
                           c('tasmax.975.mon','Warm Month Design Temperature 97.5%','degC'), ##tasmax.975.mon
                           c('tasmax.950.mon','Warm Month Design Temperature 95.0%','degC'), ##tasmax.950.mon
                           c('wb.996.mon','Warm Month Design Wet Bulb Temp 99.6%','degC'), ##wb.996.mon
                           c('wb.990.mon','Warm Month Design Wet Bulb Temp 99.0%','degC'), ##wb.990.mon
                           c('wb.975.mon','Warm Month Design Wet Bulb Temp 97.5%','degC'), ##wb.975.mon
                           c('wb.950.mon','Warm Month Design Wet Bulb Temp 95.0%','degC')) ##wb.950.mon


bc.quantile.parameters <- function(data.subsets) {

  tas.sub <- data.subsets$tas
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin
  pas.sub <- data.subsets$pas
  huss.sub <- data.subsets$huss

  tasmin.050.mon <- mean(quant.fxn(tasmin.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.050),na.rm=TRUE)
  tasmin.025.mon <- mean(quant.fxn(tasmin.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.025),na.rm=TRUE)
  tasmin.010.mon <- mean(quant.fxn(tasmin.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.010),na.rm=TRUE)
  tasmin.004.mon <- mean(quant.fxn(tasmin.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.004),na.rm=TRUE)
  tasmax.996.mon <- mean(quant.fxn(tasmax.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  tasmax.990.mon <- mean(quant.fxn(tasmax.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.990),na.rm=TRUE)
  tasmax.975.mon <- mean(quant.fxn(tasmax.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.975),na.rm=TRUE)
  tasmax.950.mon <- mean(quant.fxn(tasmax.sub$data,as.factor(format(tas.sub$dates,'%Y')),pctl=0.950),na.rm=TRUE)

  dew.point.day <- dew.point.temp(pas.sub$data,huss.sub$data)
  wet.bulb.day  <- temp.wet.bulb(tas.sub$data,dew.point.day,pas.sub$data,huss.sub$data)
  wb.996.mon <- mean(quant.fxn(wet.bulb.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)
  wb.990.mon <- mean(quant.fxn(wet.bulb.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.990),na.rm=TRUE)
  wb.975.mon <- mean(quant.fxn(wet.bulb.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.975),na.rm=TRUE)
  wb.950.mon <- mean(quant.fxn(wet.bulb.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.950),na.rm=TRUE)

  rv <- rbind(tasmin.050.mon,tasmin.025.mon,tasmin.010.mon,tasmin.004.mon,
              tasmax.996.mon,tasmax.990.mon,tasmax.975.mon,tasmax.950.mon,
              wb.996.mon,wb.990.mon,wb.975.mon,wb.950.mon)
  return(rv)
}

bc.mixed.names <- rbind(c('hdd.18.ann','Heating Degree Days (18 degC)','DD'), ##hdd.18.ann
                        c('hdd.18.3.ann','Heating Degree Days (18.3 degC)','DD'), ##hdd.18.3.ann
                        c('cdd.18.3.ann','Cooling Degree Days (18.3 degC)','DD'), ##cdd.18.3.ann
                        c('cdd.18.ann','Cooling Degree Days (18 degC)','DD'), ##cdd.18.ann
                        c('pr.ann.max','Annual Maximum Daily Precipitation','mm'), ##pr.ann.max
                        c('pr.ann.avg','Annual Average Total Precipitation','mm'),
                        c('pr_rp5','5-Year Return Period Daily Total Precipitation','mm'), ##pr.rp.year
                        c('pr_rp10','10-Year Return Period Daily Total Precipitation','mm'), ##pr.rp.year
                        c('pr_rp20','20-Year Return Period Daily Total Precipitation','mm'), ##pr.rp.year
                        c('pr_rp50','50-Year Return Period Daily Total Precipitation','mm')) ##pr.rp.year

bc.mixed.parameters <- function(data.subsets) {

  pr.sub <- data.subsets$pr 
  tas.sub <- data.subsets$tas
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin

  hdd.18.ann <- mean(hdd.18.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  hdd.18.3.ann <- mean(hdd.18.3.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  cdd.18.3.ann <- mean(cdd.18.3.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  cdd.18.ann <- mean(cdd.18.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.max <- mean(max.fxn(pr.sub$data,as.factor(format(pr.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.avg <- mean(sum.fxn(pr.sub$data,as.factor(format(pr.sub$dates,'%Y'))),na.rm=TRUE)
  pr.rp.5 <- calc.return.periods(pr.sub$data,pr.sub$dates,'pr',rp=5)  
  pr.rp.10 <- calc.return.periods(pr.sub$data,pr.sub$dates,'pr',rp=10)  
  pr.rp.20 <- calc.return.periods(pr.sub$data,pr.sub$dates,'pr',rp=20)  
  pr.rp.50 <- calc.return.periods(pr.sub$data,pr.sub$dates,'pr',rp=50)  

  rv <- c(hdd.18.ann,hdd.18.3.ann,cdd.18.3.ann,cdd.18.ann,
          pr.ann.max,pr.ann.avg,pr.rp.5,pr.rp.10,pr.rp.20,pr.rp.50)
  return(rv)
}

bc.snow.names <- rbind(c('snow_rp20','20-Year Return Period Daily Snow Load','kPa'), ##pr.rp.year
                       c('snow_rp50','50-Year Return Period Daily Snow Load','kPa')) ##pr.rp.year

bc.snow.parameters <- function(data.subsets) {
  snow.sub <- data.subsets$snow
  snow.rp20.year <- calc.return.periods(snow.sub$data*100,as.Date(snow.sub$dates),'snd',rp=20)*3.5/100
  snow.rp50.year <- calc.return.periods(snow.sub$data*100,as.Date(snow.sub$dates),'snd',rp=50)*3.5/100
  ##snow.rp.year <- snow.rp.year*3.5 ##To convert to snow load assuming snow density of 3.5kN/m^3  
  return(c(snow.rp20.year,snow.rp50.year))
}


bc.wind.names <- rbind(c('wind_rp5','5-Year Wind Pressure Return Period','Pa'), ##wind.press.5
                       c('wind_rp10','10-Year Wind Pressure Return Period','Pa'), ##wind.press.10
                       c('wind_rp50','50-Year Wind Pressure Return Period','Pa')) ##wind.press.50

bc.wind.parameters <- function(data.subsets) {

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
   title.col <- cbind(c('Building Code','Variable',data.names[,2]),c('Variable','Units',data.names[,3]))

   data.table <- rbind(data.header,range.matrix)
   data.out <- cbind(title.col,data.table)                 
   
   ts.ix <- grep('degC',data.out[,2])
   pc.ix <- grep('Percent',data.out[1,])
   data.out[ts.ix,c(pc.ix,pc.ix+1)] <- NA

   my.writedir <- '/storage/data/projects/rci/building_code/'
   write.table(data.out,file=write.file,
               sep=',',quote=FALSE,col.name=FALSE,row.name=FALSE)
}

write.site.adjusted.table <- function(means.list,tenth.list,ninetieth.list,data.names,write.file) {



   range.matrix <- cbind(means.list$past,
                         paste(tenth.list$past,ninetieth.list$past,sep=' to '),
                         means.list$site,
                         means.list$start.anom+means.list$site,
                         means.list$site*(1+means.list$start.prct/100),
                         paste(means.list$site+tenth.list$start.anom,means.list$site+ninetieth.list$start.anom,sep=' to '),
                         paste(means.list$site*(1+tenth.list$start.prct/100),means.list$site*(1+ninetieth.list$start.prct/100),sep=' to '),
                         means.list$middle.anom+means.list$site,
                         means.list$site*(1+means.list$middle.prct/100),
                         paste(means.list$site+tenth.list$middle.anom,means.list$site+ninetieth.list$middle.anom,sep=' to '),
                         paste(means.list$site*(1+tenth.list$middle.prct/100),means.list$site*(1+ninetieth.list$middle.prct/100),sep=' to '),
                         means.list$end.anom+means.list$site,
                         means.list$site*(1+means.list$end.prct/100),
                         paste(means.list$site+tenth.list$end.anom,means.list$site+ninetieth.list$end.anom,sep=' to '),
                         paste(means.list$site*(1+tenth.list$end.prct/100),means.list$site*(1+ninetieth.list$end.prct/100),sep=' to '))



    data.header  <- rbind(c('Model Past',' ',
                            'Handbook Past',
                            '2020s Site Add','Percent','Add','Percent',
                            '2050s Site Add','Percent','Add','Percent',
                            '2080s Site Add','Percent','Add','Percent'),
                        c('Average','10th%-90th%','Site',
                          'Average','10th%-90th%','Average','10th%-90th%',
                          'Average','10th%-90th%','Average','10th%-90th%',
                          'Average','10th%-90th%','Average','10th%-90th%'))
   title.col <- cbind(c('Building Code','Variable',data.names[,2]),c('Variable','Units',data.names[,3]))

   data.table <- rbind(data.header,range.matrix)
   data.out <- cbind(title.col,data.table)                 


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

my.writedir <- '/storage/data/projects/rci/building_code/granville_41st/'
##my.writedir <- '/storage/data/projects/rci/building_code/cowichan_hospital/'

##Nanaimo Hospital Coordinates
##lon.c <- -123.969357
##lat.c <- 49.184737

##Cowichan Hospital Coordinates
##lon.c <- -123.722256
##lat.c <- 48.785912

##Granville and 41st Hospital
lon.c <- -123.139690
lat.c <- 49.234305


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
        bccaq.file <- paste0(data.dir,'granville_41st_bccaq_',gcm,'_',scenario,'_pr_tx_tn.RData')
        ##bccaq.file <- paste0(data.dir,'cowichan_hospital_bccaq_',gcm,'_',scenario,'_pr_tx_tn.RData')
        if (!file.exists(bccaq.file)) {
           bccaq.data <- gather.bccaq.data(gcm,lon.c,lat.c,scenario)        
           save(bccaq.data,file=bccaq.file)
        } else {
           load(bccaq.file)
        }
        if (scenario=='rcp85') {
          snow.file <-  paste0("/storage/data/climate/downscale/RCM/CanRCM4/daily/",    
                               ##"snd_COWICHAN_HOSPITAL_77_151_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData")
                               "snd_GRANVILLE_41st_80_152_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData")
          ##paste0(data.dir,'granville_41st_bccaq_',gcm,'_',scenario,'_snow.RData')
          if (!file.exists(snow.file)) {
              browser()
             snow.data <- gather.snow.data(gcm,lon.c,lat.c,scenario)        
             save(snow.data,file=snow.file)
          } else {
             load(snow.file)
             snow.data <- rcm.data
          }
        }
        tas.data <- bccaq.data$tas
        tasmax.data <- bccaq.data$tasmax
        tasmin.data <- bccaq.data$tasmin
        pr.data <- bccaq.data$pr        
        ##snow.data <- snow.data$snow

        ##GCM other variables
        gcm.file <- paste0(data.dir,'granville_41st_gcm_',gcm,'_',scenario,'_huss_psl_uas_vas.RData')
        ##gcm.file <- paste0(data.dir,'cowichan_hospital_gcm_',gcm,'_',scenario,'_huss_psl_uas_vas.RData')
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
        past.var.subsets <- get.variables.subset(tas.data,tasmax.data,tasmin.data,
                                    psl.data,huss.data,pr.data,uas.data,vas.data,snow.data,
                                    past.int)

        start.int <- '2011-2040'
        start.var.subsets <- get.variables.subset(tas.data,tasmax.data,tasmin.data,
                                    psl.data,huss.data,pr.data,uas.data,vas.data,snow.data,
                                    start.int)
                                    
        middle.int <- '2041-2070'
        middle.var.subsets <- get.variables.subset(tas.data,tasmax.data,tasmin.data,
                                    psl.data,huss.data,pr.data,uas.data,vas.data,snow.data,
                                    middle.int)

        end.int <- '2071-2100'
        end.var.subsets <- get.variables.subset(tas.data,tasmax.data,tasmin.data,
                                    psl.data,huss.data,pr.data,uas.data,vas.data,snow.data,
                                    end.int)
                                            

        if (gcm=='HadGEM2-ES') {
          end.vars <- lapply(end.var.subsets,fix.hadley)
          end.var.subsets <- end.vars
        }

        ##------------------------------------------------------------------------------------------------
        ##------------------------------------------------------------------------------------------------
        ##ASHRAE Parameters
if (1==1) {
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

        write.combined.file <- paste(my.writedir,'ashrae.gcm.bccaq.',scenario,'.ens.combined.all.variables.granville.41st.building.code.csv',sep='')
        write.combined.table(ashrae.means,ashrae.10th,ashrae.90th,ashrae.all.names,write.combined.file)
}


        bc.means <- lapply(lapply(bc.matrix,function(x){apply(x,2,mean,na.rm=T)}),round,1)
        bc.means$start.prct <- round(bc.means$start.anom/bc.means$past*100,1)
        bc.means$middle.prct <- round(bc.means$middle.anom/bc.means$past*100,1)
        bc.means$end.prct <- round(bc.means$end.anom/bc.means$past*100,1)

        bc.10th <- lapply(lapply(bc.matrix,function(x){apply(x,2,quantile,0.1,na.rm=T)}),round,1)
        bc.10th$start.prct <- round(bc.10th$start.anom/bc.10th$past*100,1)
        bc.10th$middle.prct <- round(bc.10th$middle.anom/bc.10th$past*100,1)
        bc.10th$end.prct <- round(bc.10th$end.anom/bc.10th$past*100,1)

        bc.90th <- lapply(lapply(bc.matrix,function(x){apply(x,2,quantile,0.9,na.rm=T)}),round,1)
        bc.90th$start.prct <- round(bc.90th$start.anom/bc.90th$past*100,1)
        bc.90th$middle.prct <- round(bc.90th$middle.anom/bc.90th$past*100,1)
        bc.90th$end.prct <- round(bc.90th$end.anom/bc.90th$past*100,1)

##browser()
        write.combined.file <- paste(my.writedir,'bc.gcm.bccaq.',scenario,'.ens.combined.all.variables.granville.41st.building.code.csv',sep='')
        ##write.combined.table(bc.means,bc.10th,bc.90th,bc.all.names,write.combined.file)
        ##write.combined.file <- paste(my.writedir,'bc.gcm.bccaq.',scenario,'.ens.formatted.all.variables.cowichan.hospital.building.code.csv',sep='')
        write.combined.table(bc.means,bc.10th,bc.90th,bc.all.names,write.combined.file)
        write.adjusted.file <- paste(my.writedir,'bc.gcm.bccaq.',scenario,'.ens.formatted.site.adjusted.granville.41st.building.code.csv',sep='') 
        write.site.adjusted.table(bc.means,bc.10th,bc.90th,bc.all.names,write.adjusted.file)

}



