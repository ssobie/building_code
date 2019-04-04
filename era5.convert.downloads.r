##Script to convert the downloaded ERA5 data into useful
##file formats and structures

library(ncdf4)
library(udunits2)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R')

##-------------------------------------------------------------------------------------
##Uses Gregorian Calendar 
era5_time_series <- function(freq,dates) {

  calendar <- 'gregorian'
  origin <- '1950-01-01 00:00:00'
  pcict.origin <- as.PCICt(origin,cal=calendar)
  time.values <- dates-pcict.origin  

  ##Daily 
  daily.series <- format(round.PCICt(dates,freq),'%Y-%m-%d')
  daily.values <- as.numeric(time.values/86400 ) 
  daily.time <- list(calendar=calendar,
                    freq=freq,
                    units=paste0('days since ',origin),
                    long_name='time',
                    standard_name='time',
                    values=daily.values,
                    series=daily.series)

  ##Hourly 
  hourly.series <- format(round.PCICt(dates,freq),'%Y-%m-%d %H:%M:%S')
  hourly.values <- as.numeric(time.values/3600 ) 
  hourly.time <- list(calendar=calendar,
                    freq=freq,
                    units=paste0('hours since ',origin),
                    long_name='time',
                    standard_name='time',
                    values=hourly.values,
                    series=hourly.series)

  rv <- switch(freq,
               day=daily.time,
               hour=hourly.time)
  return(rv) 
}

##-------------------------------------------------------------------------------------
##Global ERA5 Attributes
get_global_atts <- function(freq) {
  global.atts <- list(institution="European Centre for Medium-Range Weather Forecasts",
                   contact="ECMWF",
                   Conventions="CF-1.6",
                   institute_id ="ECMWF",
                   domain='British Columbia',
                   creation_date=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%Z'),
                   frequency=freq,
                   product="reanalysis",
                   modeling_realm="atmos",
                   project_id='ERA5',
                   references="Copernicus Climate Change Service (C3S) (2017): ERA5: Fifth generation of ECMWF atmospheric reanalyses of the global climate .                         Copernicus Climate Change Service Climate Data Store (CDS), 28 March 2019. https://cds.climate.copernicus.eu/cdsapp#!/home")

}
##-------------------------------------------------------------------------------------
get_variable_units <- function(var.name) {

  units <- list(prhour="kg m-2 h-1",prday="kg m-2 day-1",
                tasmax='degC',tasmin='degC',tasday='degC',tashour='degC',tasrange='degC',
                tasskew='')             
  return(units[[var.name]])
}
##-------------------------------------------------------------------------------------

get_variable_specific_atts <- function(var.name) {

  pr.hour.atts <- list(units = get_variable_units('prhour'))
  pr.day.atts <- list(units = get_variable_units('prday'))
  tasmax.atts <- list(long_name = "Daily Maximum Near-Surface Air Temperature",
                      cell_methods = "time: maximum",units=get_variable_units('tasmax'))
  tasmin.atts <- list(long_name = "Daily Minimum Near-Surface Air Temperature",
                      cell_methods = "time: minimum",units=get_variable_units('tasmin'))
  tasday.atts <- list(long_name = "Daily Average Near-Surface Air Temperature",
                      cell_methods = "time: mean",units=get_variable_units('tasday'))
  tashour.atts <- list(long_name = "Hourly Average Near-Surface Air Temperature",
                       cell_methods = "time: mean",units=get_variable_units('tashour'))
  tasrange.atts <- list(long_name = "Daily Near-Surface Air Temperature Range",
                       cell_methods = "time: range",units=get_variable_units('tasrange'))
  tasskew.atts <- list(long_name = "Daily Near-Surface Air Temperature Skewness",
                       cell_methods = "time: skewness",units=get_variable_units('tasskew'))

  var.atts <- switch(var.name,
                     prhour=pr.hour.atts,
                     prday=pr.day.atts,
                     tasmax=tasmax.atts,
                     tasmin=tasmin.atts,
                     tasday=tasday.atts,
                     tashour=tashour.atts,
                     tasrange=tasrange.atts,
                     tasskew=tasskew.atts)
  rv <- list(var=var.atts)
  return(rv)
}

##-------------------------------------------------------------------------------------
get_standard_atts <- function(var.name) {
  lon.atts <- list(standard_name="longitude",long_name = "longitude",
                   units = "degrees_east",axis = "X")

  lat.atts <- list(standard_name="latitude",long_name = "latitude",
                   units = "degrees_north",axis = "Y")

  pr.atts <- list(standard_name = "total_precipitation",
                  long_name = "Precipitation",
                  missing_value = 1.e+20,
                  cell_methods = "time: sum")

  tas.atts <- list(standard_name = "air_temperature",
                   missing_value = 1.e+20)

  var.atts <- switch(var.name,
                     pr=pr.atts,
                     tas=tas.atts)

  rv <- list(lon=lon.atts,
             lat=lat.atts,
             var=var.atts)
  return(rv)
}

add_attributes_ncdf <- function(var.info, time, nc) {

  standard.atts <- get_standard_atts(var.info$type)
  variable.atts <- get_variable_specific_atts(var.info$name)
  print('Lon names')
  lon.names <- names(standard.atts$lon)
  for (j in 1:length(standard.atts$lon))
    ncatt_put(nc,varid='lon',attname=lon.names[j],attval=standard.atts$lon[[j]])
  print('Lat names')
  lat.names <- names(standard.atts$lat)
  for (j in 1:length(standard.atts$lat))
    ncatt_put(nc,varid='lat',attname=lat.names[j],attval=standard.atts$lat[[j]])

  print('Standard names')
  var.names <- names(standard.atts$var)
  for (j in 1:length(standard.atts$var))
    ncatt_put(nc,varid=var.info$filevar,attname=var.names[j],attval=standard.atts$var[[j]])
  print('Variable names')
  var.names <- names(variable.atts$var)
  for (j in 1:length(variable.atts$var))
    ncatt_put(nc,varid=var.info$filevar,attname=var.names[j],attval=variable.atts$var[[j]])
 
  print('Time atts')
  ##Time attributes
  ncatt_put(nc,varid='time',attname='units',attval=time$units)
  ncatt_put(nc,varid='time',attname='long_name',attval=time$long_name)
  ncatt_put(nc,varid='time',attname='standard_name',attval=time$standard_name)
  ncatt_put(nc,varid='time',attname='calendar',attval=time$calendar)

  print('Global atts')
  ##Global Attributes
  global.atts <- get_global_atts(freq=time$freq)
  global.names <- names(global.atts)
  for (g in 1:length(global.atts))
    ncatt_put(nc,varid=0,attname=global.names[g],attval=global.atts[[g]])

  ##Clear extraneous history
  ncatt_put(nc,varid=0,attname='history',attval='')

}


##-------------------------------------------------------------------------------------
make_era5_netcdf <- function(var.info,base.file,dir) {

  nc <- nc_open(paste0(dir,base.file),write=FALSE)
  lon <- ncvar_get(nc,'longitude')
  lat <- ncvar_get(nc,'latitude')  

  time <- era5_time_series(var.info$freq,var.info$dates)

  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time$units, time$values,
                      unlim=FALSE, calendar=time$calendar)

  var.geog <- ncvar_def(var.info$filevar, units=get_variable_units(var.info$name),
                        dim=list(x.geog, y.geog,t.geog),
                        missval=1.e+20)
  new.nc <- nc_create(paste(dir,var.info$file,sep=''),var.geog)

  add_attributes_ncdf(var.info, time, new.nc)                        
  ncvar_put(new.nc,'lon',lon)
  ##Invert Lat values
  n.lat <- length(lat)
  ncvar_put(new.nc,'lat',lat[(n.lat:1)])
  
  nc_close(nc)
  nc_close(new.nc)
  
}

##-------------------------------------------------------------------------------------

make_daily_series <- function(year,input.data,agg.fxn) {

   year.dates <- seq(from=as.PCICt(paste0(year,'-01-01 00:00:00'),cal=calendar),
                     by='hour',
                     to=as.PCICt(paste0(year,'-12-31 23:00:00'),cal=calendar))

   day.fac <- as.factor(format(year.dates,'%Y-%m-%d'))
   day.time <- as.Date(levels(day.fac))
  
   input.agg <- aperm(apply(input.data,c(1,2),function(x,y){tapply(x,y,agg.fxn)},day.fac),c(2,3,1))

   return(input.agg)
}

 
##-------------------------------------------------------------------------------------
concatenate_series <- function(var.name,year,vnc,var.info,input.data) {
   all.dates <- var.info$dates
   freq <- var.info$freq
   year.dates <- seq(from=as.PCICt(paste0(year,'-01-01 00:00:00'),cal=calendar),
                     by=freq,
                     to=as.PCICt(paste0(year,'-12-31 23:00:00'),cal=calendar))
   yix <- which(all.dates %in% year.dates)
   yst <- head(yix,1)
   yen <- tail(yix,1)
   yct <- yen-yst+1

   ncvar_put(vnc,var.info$filevar,input.data,start=c(1,1,yst),count=c(-1,-1,yct))
}
 
##-------------------------------------------------------------------------------------


##*************************************************************************************

read.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/ERA5/'
tmp.dir <- '/local_temp/ssobie/era5/'
if (!file.exists(tmp.dir))
   dir.create(tmp.dir,recursive=TRUE)

##-------------------------------------------------------------------------------------
##From Hourly Temperature

##File for metadata attributes
base.file <- '2m_temperature_hour_ERA5_BC_1995.nc'
file.copy(from=paste0(read.dir,base.file),to=tmp.dir,overwrite=TRUE)

era5.var <- 't2m'
calendar <- 'gregorian'
hour.dates <- seq(from=as.PCICt('1980-01-01 00:00:00',cal=calendar),by='hour',to=as.PCICt('2018-12-31 23:00:00',cal=calendar))
day.dates <- seq(from=as.PCICt('1980-01-01',cal=calendar),by='day',to=as.PCICt('2018-12-31',cal=calendar))
years <- 1980:2018

##day.time <- era5_time_series('day',day.dates)
##hour.time <- era5_time_series('hour',hour.dates)


##Create netcdf write files
var.list <- c('tashour','tasday','tasmax','tasmin','tasrange','tasskew')
info.list <- list(tashour=list(name='tashour',filevar='tas',type='tas',freq='hour',dates=hour.dates,file='tas_hour_ERA5_BC_19800101-20181231.nc'),  
                  tasday=list(name='tasday',filevar='tas',type='tas',freq='day',dates=day.dates,file='tas_day_ERA5_BC_19800101-20181231.nc'),  
                  tasmax=list(name='tasmax',filevar='tasmax',type='tas',freq='day',dates=day.dates,file='tasmax_day_ERA5_BC_19800101-20181231.nc'),                
                  tasmin=list(name='tasmin',filevar='tasmin',type='tas',freq='day',dates=day.dates,file='tasmin_day_ERA5_BC_19800101-20181231.nc'),                
                  tasrange=list(name='tasrange',filevar='tasrange',type='tas',freq='day',dates=day.dates,file='tasrange_day_ERA5_BC_19800101-20181231.nc'),
                  tasskew=list(name='tasskew',filevar='tasskew',type='tas',freq='day',dates=day.dates,file='tasskew_day_ERA5_BC_19800101-20181231.nc'))


##----------------------------
##Create derived files
ncs.list <- vector(mode='list',length=length(var.list))
names(ncs.list) <- var.list
for (var.name in var.list) {
   print(paste0('Making file for ',var.name))
   var.info <- info.list[[var.name]]
##   make_era5_netcdf(var.info,base.file=base.file,
##                    dir=tmp.dir)
   ncs.list[[var.name]] <- nc_open(paste0(tmp.dir,var.info$file),write=TRUE)
 
}

##---------------------------

for (year in years) {
  print(year)
  year.file <- paste0('2m_temperature_hour_ERA5_BC_',year,'.nc')
  file.copy(from=paste0(read.dir,year.file),to=tmp.dir,overwrite=TRUE)

  ync <- nc_open(paste0(tmp.dir,year.file))
  tas.raw <- ncvar_get(ync,'t2m')
  tas.inv <- ud.convert(tas.raw,'K','degC')

  ##Invert tas to fix latitude order
  n.lat <- ync$dim$latitude$len
  tas <- tas.inv[,(n.lat:1),]

  ##Concatenate Hourly temperatures
  ##concatenate_series('tashour',year,ncs.list[['tashour']],info.list[['tashour']],tas)

  ##Calculate Daily Maximum Temperature
  tasmax <- make_daily_series(year,tas,max)
  concatenate_series('tasmax',year,ncs.list[['tasmax']],info.list[['tasmax']],tasmax)
  tasmin <- make_daily_series(year,tas,min)
  concatenate_series('tasmin',year,ncs.list[['tasmin']],info.list[['tasmin']],tasmin)
  tas.day <- make_daily_series(year,tas,mean)
  concatenate_series('tasday',year,ncs.list[['tasday']],info.list[['tasday']],tas.day)
  tasrange <- tasmax-tasmin
  concatenate_series('tasrange',year,ncs.list[['tasrange']],info.list[['tasrange']],tasrange)
  tasskew <- (tas.day - tasmin)/tasrange
  concatenate_series('tasskew',year,ncs.list[['tasskew']],info.list[['tasskew']],tasskew)

  nc_close(ync)
  file.remove(paste0(tmp.dir,year.file))

}


for (var.name in var.list) {
##  nc_close(ncs.list[[var.name]])
  print(info.list[[var.name]][['file']])
  file.copy(from=paste0(tmp.dir,info.list[[var.name]][['file']]),to=paste0(read.dir,'concat/'),overwrite=TRUE)
}





##Calculate TASMAX
##Calculate TASMIN
##Calculate TAS
##Calculate TAS RANGE
##Calculate TAS SKEW
##Concatenate daily values into single file



