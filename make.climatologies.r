##Script to add attributes to the aggregated downscaled output.

library(ncdf4)
library(PCICt)

sub.by.time <- function(input.file,output.file,interval,read.dir,write.dir) {

  nc <- nc_open(paste(read.dir,input.file,sep=''))
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*86400
  years <- format(time.series,'%Y')
  yrs <- strsplit(interval,'-')[[1]]
  
  st <- head(grep(yrs[1],years),1)-1
  en <- tail(grep(yrs[2],years),1)-1
  sub.file <- output.file
  system(paste('ncks --overwrite -d time,',st,',',en,' ',read.dir,input.file,' ',write.dir,sub.file,sep=''))
  nc_close(nc)
  rv <- paste0(write.dir,sub.file)
  return(rv)
}


run.annual.climatologies <- function(gcm,scenario) {

  intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')

  proj.dir <-  '/storage/data/climate/downscale/CMIP5/building_code/'

  print(gcm)
  read.dir <- paste(proj.dir,gcm,'/',sep='')
  write.dir  <- paste0(read.dir,'climatologies/')
  if (!file.exists(write.dir)) {
     dir.create(write.dir,recursive=T)
  }

  ##------------------------------------------------------------------------------
  ##Maximum Temperature
  all.files <- list.files(path=read.dir,pattern='wetbulb_day')
  day.file <- all.files[grep(scenario,all.files)]
  
  int.file <- gsub(pattern='wetbulb_day',replacement='wetbulb_interval',day.file)  
  all.file <- gsub(pattern='wetbulb_day',replacement='wetbulb_all_quantile_010',day.file)  
  ann.file <- gsub(pattern='wetbulb_day',replacement='wetbulb_annual_quantile_010',day.file)  
  ##jul.file <- gsub(pattern='wetbulb_day',replacement='wetbulb_month_010',day.file)  
  ##jan.file <- gsub(pattern='wetbulb_day',replacement='wetbulb_month_001',day.file)  
  ##avg.file <- gsub(pattern='wetbulb_annual',replacement='wetbulb_average_annual_climatology',ann.file)  
  for (interval in intervals) {     
     print(paste0('dp ',interval))   
     sub.file <- sub.by.time(day.file,int.file,interval=interval,read.dir,write.dir)
     system(paste('cdo -s -O yearpctl,0.010 ',sub.file,' -yearmin ',sub.file,' -yearmax ',sub.file,' ',
                   write.dir,gsub(pattern='19500101-21001231',replacement=interval,all.file),sep=''))
     ##system(paste('cdo -s -O ymonpctl,0.010 ',sub.file,' -ymonmin ',sub.file,' -ymonmax ',sub.file,' ',
     ##              write.dir,gsub(pattern='19500101-21001231',replacement=interval,jul.file),sep=''))
     ##system(paste('cdo -s -O ymonpctl,0.01 ',sub.file,' -ymonmin ',sub.file,' -ymonmax ',sub.file,' ',
     ##              write.dir,gsub(pattern='19500101-21001231',replacement=interval,jan.file),sep=''))
     system(paste('cdo -s -O timmean ',write.dir,gsub(pattern='19500101-21001231',replacement=interval,all.file),' ',write.dir,gsub(pattern='19500101-21001231',replacement=interval,ann.file),sep=''))

     file.remove(sub.file)         
     ##file.remove(paste0(write.dir,int.file))         
     file.remove(paste0(write.dir,gsub(pattern='19500101-21001231',replacement=interval,all.file)))
  }

}


  

##************************************************************************

gcm.list <- c('ACCESS1-0',
              'CanESM2',
              'CNRM-CM5',
              'CSIRO-Mk3-6-0',
              'GFDL-ESM2G',
              'HadGEM2-CC',
              'HadGEM2-ES',
              'inmcm4',
              'MIROC5',
              'MRI-CGCM3')
scenario <- 'rcp85'

##Annual
for (gcm in gcm.list) {
  run.annual.climatologies(gcm,scenario)


}
