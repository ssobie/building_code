##Contains functions to compute Building Code values from the BC
##Building Code and the American ASHRAE Handbook

library(ncdf4)
library(PCICt)

source('/storage/home/ssobie/code/repos/building_code/building.code.fxns.r')
source('/storage/home/ssobie/code/repos/building_code/convert.ashrae.to.metric.r')

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

  snd.mon <- get.data.subset(snd.data,dates,interval) 

  dates.mon <- get.data.subset(dates,dates,interval)

  rv <- list(tas=tas.mon,tasmax=tasmax.mon,tasmin=tasmin.mon,
             pas=pas.mon,huss=huss.mon,pr=pr.mon,
             uas=uas.mon,vas=vas.mon,snd=snd.mon,dates=dates.mon)
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
  
  rv <- c(wspd.high,mcdb.wspd.high,wspd.mcdb.low,db.pcwd.low,wspd.mcdb.high,db.pcwd.high,wspd.2.5.ann,wspd.1.ann)
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
                               c('Cooling Dry Bulb Temperature 99.6%','degC'),  ##db.cool.ann
                               c('Average Annual Max Dry Bulb Temperature','degC'),  ##db.max.ann
                               c('Average Annual Min Dry Bulb Temperature','degC'))  ##db.min.ann

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


ashrae.humidity.names <- rbind(c('Cold Dew Point Temperature 0.4% degC','degC'),            ##dew.point.low
                               c('Cold Humidity Ratio 0.4%','Ratio'),                       ##hum.ratio.low
                               c('Mean Coincident Dry Bulb (Dew Point 0.4%) degC','degC'),  ##dew.mcdb.low
                               c('Mean Coincident Wet Bulb (Dry Bulb 99.6%) degC','degC'),  ##db.mcwb.ann
                               c('Evaporation Wet Bulb 99.6% degC','degC'),                 ##wb.evap.ann
                               c('Mean Coincident Dry Bulb (Wet Bulb 99.6%) degC','degC'),  ##wb.mcdb.ann
                               c('Warm Dew Point Temperature 99.6% degC','degC'),           ##dew.point.high
                               c('Warm Humidity Ratio 99.6%','Ratio'),                      ##hum.ratio.high
                               c('Mean Coincident Dry Bulb (Dew Point 99.6%) degC','degC'), ##dew.mcdb.high
                               c('Enthalpy 99.6% kJ/kg','kJ/kg'),                           ##enth.ann
                               c('Mean Coincident Dry Bulb (Enthalpy 99.6%) degC','degC'),  ##enth.mcdb.high
                               c('Annual Maximum Wet Bulb degC','degC'))                    ##wb.max.ann       


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
  enth.ann <- mean(quant.fxn(enth.day,as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE)/1000 ##to convert to kJ
  enth.mcdb.high <- mean(conditional.mon.fxn(main=tas.sub$data,second=enth.day,fac=as.factor(format(tas.sub$dates,'%Y')),pctl=0.996),na.rm=TRUE) 
  
  wb.max.ann <- mean(max.fxn(wet.bulb.day,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  browser()
  rv <- c(dew.point.low,hum.ratio.low,dew.mcdb.low,db.mcwb.ann,wb.evap.ann,wb.mcdb.ann,
          dew.point.high,hum.ratio.high,dew.mcdb.high,enth.ann,enth.mcdb.high,wb.max.ann)
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

bc.mixed.names <- rbind(c('Heating Degree Days 18 degC','DD'), ##hdd.18.ann
                        c('Annual Maximum Daily Precipitation','mm'), ##pr.ann.max
                        c('Annual Average Total Precipitation','mm'),
                        c('50-Year Return Period Daily Total Precipitation','mm')) ##pr.rp.year


bc.mixed.parameters <- function(data.subsets) {

  pr.sub <- data.subsets$pr 
  tas.sub <- data.subsets$tas
  tasmax.sub <- data.subsets$tasmax
  tasmin.sub <- data.subsets$tasmin


  hdd.18.ann <- mean(hdd.18.fxn(tas.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.max <- mean(max.fxn(pr.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  pr.ann.avg <- mean(sum.fxn(pr.sub$data,as.factor(format(tas.sub$dates,'%Y'))),na.rm=TRUE)
  pr.rp.year <- calc.return.periods(pr.sub$data,tas.sub$dates,'pr',rp=50)  

  rv <- c(hdd.18.ann,pr.ann.max,pr.ann.avg,pr.rp.year)
  return(rv)
}

bc.snow.names <- rbind(c('50-Year Return Period Daily Snow Load','kPa'), ##pr.rp.year
                       c('50-Year Return Period Daily Snow Load','kPa')) ##pr.rp.year
bc.snow.parameters <- function(data.subsets) {
  snow.sub <- data.subsets$snd
  snow.rp.year <- calc.return.periods(snow.sub$data*100,snow.sub$dates,'snd',rp=50)
  snow.rp.year <- snow.rp.year/35 ##To convert to snow load assuming snow density of 3.5kN/m^3  
  return(snow.rp.year)
}


bc.wind.names <- rbind(c('5-Year Wind Pressure Return Period','kPa'), ##wind.press.5
                       c('10-Year Wind Pressure Return Period','kPa'), ##wind.press.10
                       c('50-Year Wind Pressure Return Period','kPa')) ##wind.press.50

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

  rv <- c(wind.press.5,wind.press.10,wind.press.50)
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
                 paste('Site Past (',var.units,')',sep=''),
                 paste('Past (',var.units,')',sep=''),
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


##************************************************************************
##************************************************************************
##Load Data

##RCM Only Data
##source('/storage/home/ssobie/code/repos/building_code/rcm.build.code.data.r',chdir=TRUE)
#source('/storage/home/ssobie/code/repos/building_code/gcm.bccaq.build.code.data.r',chdir=TRUE)

read.dir <- '/storage/data/climate/downscale/BCCAQ2/CanRCM4/'
cell <- '76_150'
#load(paste(read.dir,'snd_NANAIMO_',cell,'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
#snd.data <- data


##-----------------------------------------------------------------------------
##Monthly Parameters
##Past
past.int <- '1971-2000'
past.var.subsets <- get.variables.subset(tas.data$data,tasmax.data$data,tasmin.data$data,
                                    psl.data$data,huss.data$data,pr.data$data,uas.data$data,vas.data$data,snd.data,
                                    tasmax.data$time,past.int)

start.int <- '2011-2040'
start.var.subsets <- get.variables.subset(tas.data$data,tasmax.data$data,tasmin.data$data,
                                    psl.data$data,huss.data$data,pr.data$data,uas.data$data,vas.data$data,snd.data,
                                    tasmax.data$time,start.int)

middle.int <- '2041-2070'
middle.var.subsets <- get.variables.subset(tas.data$data,tasmax.data$data,tasmin.data$data,
                                    psl.data$data,huss.data$data,pr.data$data,uas.data$data,vas.data$data,snd.data,
                                    tasmax.data$time,middle.int)

end.int <- '2071-2100'
end.var.subsets <- get.variables.subset(tas.data$data,tasmax.data$data,tasmin.data$data,
                                    psl.data$data,huss.data$data,pr.data$data,uas.data$data,vas.data$data,snd.data,
                                    tasmax.data$time,end.int)

##------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------

##Nanaimo Variables
ashrae.vars <-  list(c(11.2,'degF'), ##Cold Dewpoint Temp 0.4% degF
                     c(9.8,'ratio'),  ##Cold Humidity Ratio 0.4% gr/lb
                     c(28.4,'degF'), ##MCDB Dewpoint 0.4% degF
                     c(63.9,'degF'), ##MCWB Dry Bulb 99.6% degF
                     c(65.2,'degF'), ##Evap Wet Bulb 99.6% degF
                     c(76.8,'degF'), ##MCDB Wet Bulb 99.6% degF
                     c(60.8,'degF'), ##Warm Dewpoint 99.6% degF (Dehumidification)
                     c(79.9,'ratio'), ##Warm Humidity Ratio 99.6% gr/lb
                     c(68.9,'degF'),  ##MCDB Depoint 99.6% degF
                     c(30.1,'btu.lb'), ##Enthalpy 99.6% Btu/lb
                     c(77.0,'degF'),  ##MCDB Enthalpy 99.6% degF
                     c(75.2,'degF'), ##Annual Max Wet Bulb degF
                     c(23.4,'degF'), ##Heating Dry Bulb 0.4% degF
                     c(80.3,'degF'), ##Cooling Dry Bulb 99.6% degF
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
                     c(25.9,'mph'),  ##Extreme Annual Wind Speed 97.5% mph
                     c(30.2,'mph'),  ##Extreme Annual Wind Speed 99% mph
                     c(89.9,'degF'), ##5-Year Return Period Max Temp degF
                     c(17.0,'degF'), ##5-Year Return Period Min Temp degF
                     c(NA,'inch')) ##5-Year Return Period Daily Precip inches

ashrae.si.versions <- unlist(lapply(ashrae.vars,convert.to.si))

my.writedir <- '/storage/data/projects/rci/building_code/'

##ASHRAE Parameters
ashrae.humidity <- get.build.code.anomalies(ashrae.humidity.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
##write.file <- paste(my.writedir,'ashrae.humidity.variables.nanaimo.building.code.csv',sep='')
##write.annual.table(ashrae.humidity,ashrae.humidity.names,write.file)

ashrae.dry.temp <- get.build.code.anomalies(ashrae.annual.dry.temp.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
##write.file <- paste(my.writedir,'ashrae.dry.temp.variables.nanaimo.building.code.csv',sep='')
##write.annual.table(ashrae.dry.temp,ashrae.dry.temp.names,write.file)

ashrae.precip   <- get.build.code.anomalies(ashrae.annual.precip.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
##write.file <- paste(my.writedir,'ashrae.precip.variables.nanaimo.building.code.csv',sep='')
##write.annual.table(ashrae.precip,ashrae.precip.names,write.file)

ashrae.wind     <- get.build.code.anomalies(ashrae.annual.wind.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
##write.file <- paste(my.writedir,'ashrae.wind.variables.nanaimo.building.code.csv',sep='')
##write.annual.table(ashrae.wind,ashrae.wind.names,write.file)

ashrae.rp       <- get.build.code.anomalies(ashrae.return.periods.5,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
##write.file <- paste(my.writedir,'ashrae.return.period.variables.nanaimo.building.code.csv',sep='')
##write.annual.table(ashrae.rp,ashrae.rp.names,write.file)

ashrae.all.parameters <- mapply(c,ashrae.humidity,ashrae.dry.temp,ashrae.precip,ashrae.wind,ashrae.rp,SIMPLIFY=FALSE)
ashrae.all.names <- rbind(ashrae.humidity.names,
                          ashrae.dry.temp.names,
                          ashrae.precip.names,
                          ashrae.wind.names,
                          ashrae.rp.names)

ashrae.all.parameters <- c(list(site=round(ashrae.si.versions,1)),ashrae.all.parameters)

write.file <- paste(my.writedir,'ashrae.all.variables.nanaimo.building.code.gb.csv',sep='')
write.annual.table(ashrae.all.parameters,ashrae.all.names,write.file)


##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
##BC Parameters

##Nanaimo Variables
bc.si.versions <-  c(-6.0, ##'degC'), ##Cold Month Design Temperature 2.5% degC
                     -8.0, ##'degC'),  ##Cold Month Design Temperature 1.0% degC
                     27,   ##'degC'), ##Warm Month Design Temperature 97.5% degC
                     19,   ##'degC'), ##Warm Month Wet Bulb Temperature 97.5% degC
                     3000, ##'degC'), ##Degree Days below 18 degC
                     NA,   ##'mm'), ##Annual Max Daily Precipitation mm
                     1050, ##'mm'), ##Annual Average Total Precipitation mm
                     91,   ##'mm'), ##50-Year return period Daily Precipitation
                     2.3,  ##'kPa'),  ##50-Year return period Snow Load
                     0.4,  ##'kPa'),  ##50-Year return period Rain Load
                     0.39, ##'kPa'),  ##5-Year return period Wind Load
                     0.50) ##'kPa'))  ##10-Year return period Wind Load

bc.quantiles <- get.build.code.anomalies(bc.quantile.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
##write.file <- paste(my.writedir,'bc.quantile.variables.nanaimo.building.code.csv',sep='')
##write.annual.table(bc.quantiles,bc.quantile.names,write.file)

bc.mixed <- get.build.code.anomalies(bc.mixed.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
##write.file <- paste(my.writedir,'bc.mixed.variables.nanaimo.building.code.csv',sep='')
##write.annual.table(bc.mixed,bc.mixed.names,write.file)

bc.snow <- get.build.code.anomalies(bc.snow.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,middle.var.subsets)
##write.file <- paste(my.writedir,'bc.snow.variables.nanaimo.building.code.csv',sep='')
##write.annual.table(bc.snow,bc.snow.names,write.file)


bc.wind <- get.build.code.anomalies(bc.wind.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
##write.file <- paste(my.writedir,'bc.wind.variables.nanaimo.building.code.csv',sep='')
##write.annual.table(bc.wind,bc.wind.names,write.file)
                      
bc.all.parameters <- mapply(c,bc.quantiles,bc.mixed,bc.snow,bc.wind,SIMPLIFY=FALSE)
bc.all.names <- rbind(bc.quantile.names,
                          bc.mixed.names,
                          bc.snow.names[1,],
                          bc.wind.names)

bc.all.parameters <- c(list(site=round(bc.si.versions,2)),bc.all.parameters)

write.file <- paste(my.writedir,'bc.all.variables.nanaimo.building.code.gb.csv',sep='')
write.annual.table(bc.all.parameters,bc.all.names,write.file)

      



