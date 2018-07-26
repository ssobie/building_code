##Contains functions to compute Building Code values from the BC
##Building Code and the American ASHRAE Handbook

library(ncdf4)
library(PCICt)

source('/storage/home/ssobie/code/repos/building_code/bc.building.code.fxns.r')
source('/storage/home/ssobie/code/repos/building_code/gcm.bccaq.build.code.data.r',chdir=TRUE)

##------------------------------------------------------------------
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

        
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
        ##BC Parameters
        bc.quantiles <- get.build.code.anomalies(bc.quantile.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
        bc.mixed <- get.build.code.anomalies(bc.mixed.parameters,past.var.subsets,start.var.subsets,middle.var.subsets,end.var.subsets)
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

        write.combined.file <- paste(my.writedir,'TEST.bc.gcm.bccaq.',scenario,'.ens.combined.all.variables.nanaimo.building.code.csv',sep='')
        write.combined.table(bc.means,bc.10th,bc.90th,bc.all.names,write.combined.file)


}



