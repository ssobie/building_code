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

  units <- list(pr="kg m-2 day-1",
                tasmax='degC',tasmin='degC',tasday='degC',tashour='degC',tasrange='degC',
                tasskew='',sp='Pa',uwind='m s-1',vwind='m s-1',wspd='m s-1',dewpoint='degC',
                tcc='fraction',lcc='fraction',rain='kg m-2',
                ssrd='kJ m-2',strd='kJ m-2',tisr='kJ m-2',fdir='kJ m-2')
  return(units[[var.name]])
}
##-------------------------------------------------------------------------------------

get_variable_specific_atts <- function(var.name) {

  pr.atts <- list(units = get_variable_units('pr'))
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
  sp.atts <- list(long_name = "Surface pressure",
                       cell_methods = "time: skewness",units=get_variable_units('sp'))
  uwind.atts <- list(long_name = "10 metre U wind component",
                     units=get_variable_units('uwind'))
  vwind.atts <- list(long_name = "10 metre V wind component",
                     units=get_variable_units('vwind'))
  wspd.atts <- list(long_name = "10 metre wind speed",
                     units=get_variable_units('wspd'))
  dewpoint.atts <- list(long_name = "2 metre dewpoint temperature",
                        units=get_variable_units('dewpoint'))
  totalcloud.atts <- list(long_name = "Total cloud cover",
                        units=get_variable_units('totalcloud'))
  lowcloud.atts <- list(long_name = "Low cloud cover",
                        units=get_variable_units('lowcloud'))
  totalrain.atts <- list(long_name = "Total Column Rain Water",
                        units=get_variable_units('totalrain'))


  var.atts <- switch(var.name,
                     pr=pr.atts,
                     tasmax=tasmax.atts,
                     tasmin=tasmin.atts,
                     tasday=tasday.atts,
                     tashour=tashour.atts,
                     tasrange=tasrange.atts,
                     tasskew=tasskew.atts,
                     sp=sp.atts,
                     uwind=uwind.atts,
                     vwind=vwind.atts,
                     wspd=wspd.atts,
                     dewpoint=dewpoint.atts,
                     tcc=totalcloud.atts,
                     lcc=lowcloud.atts,
                     rain=totalrain.atts)
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

  sp.atts <- list(standard_name = "surface_pressure",
                   missing_value = 1.e+20)

  uwind.atts <- list(standard_name = "10 metre U wind component",
                   missing_value = 1.e+20)
  vwind.atts <- list(standard_name = "10 metre U wind component",
                   missing_value = 1.e+20)
  wspd.atts <- list(standard_name = "10 metre wind speed",
                   missing_value = 1.e+20)

  dewpoint.atts <- list(standard_name = "2 metre dewpoint temperature",
                   missing_value = 1.e+20)
  totalcloud.atts <- list(standard_name = "Cloud cover fraction",
                   missing_value = 1.e+20)
  lowcloud.atts <- list(standard_name = "Cloud cover fraction",
                   missing_value = 1.e+20)
  totalrain.atts <- list(standard_name = "Total Column Rain Water",
                   missing_value = 1.e+20)

  var.atts <- switch(var.name,
                     pr=pr.atts,
                     tas=tas.atts,
                     sp=sp.atts,
                     uwind=uwind.atts,
                     vwind=vwind.atts,
                     wspd=wspd.atts,
                     dewpoint=dewpoint.atts,
                     totalcloud=totalcloud.atts,
                     lowcloud=lowcloud.atts,
                     totalrain=totalrain.atts)

  rv <- list(lon=lon.atts,
             lat=lat.atts,
             var=var.atts)
  return(rv)
}

add_attributes_ncdf <- function(var.name,var.type, time, nc) {

  standard.atts <- get_standard_atts(var.type)
  variable.atts <- get_variable_specific_atts(var.name)
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
    ncatt_put(nc,varid=var.name,attname=var.names[j],attval=standard.atts$var[[j]])
  print('Variable names')
  var.names <- names(variable.atts$var)
  for (j in 1:length(variable.atts$var))
    ncatt_put(nc,varid=var.name,attname=var.names[j],attval=variable.atts$var[[j]])
 
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
make_daily_era5_netcdf <- function(var.name,var.type,hour.dates,base.file,dir) {

  daily.file <- gsub('hour','day',base.file)

  nc <- nc_open(paste0(dir,base.file),write=FALSE)
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')  

  time <- era5_time_series('day',hour.dates)

  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time$units, time$values,
                      unlim=FALSE, calendar=time$calendar)

  var.geog <- ncvar_def(var.name, units=get_variable_units(var.name),
                        dim=list(x.geog, y.geog,t.geog),
                        missval=1.e+20)
  new.nc <- nc_create(paste(dir,daily.file,sep=''),var.geog)

  add_attributes_ncdf(var.name,var.type, time, new.nc)                        
  ncvar_put(new.nc,'lon',lon)
  ncvar_put(new.nc,'lat',lat)
  
  nc_close(nc)
  nc_close(new.nc)
  return(daily.file)  

}

##-------------------------------------------------------------------------------------

make_daily_series <- function(year.dates,input.data,agg.fxn) {
   day.fac <- as.factor(format(year.dates,'%Y-%m-%d'))
   day.time <- as.Date(levels(day.fac)) 
   input.agg <- aperm(apply(input.data,c(1,2),function(x,y){tapply(x,y,agg.fxn)},day.fac),c(2,3,1))
   return(input.agg)
}
 
##-------------------------------------------------------------------------------------
##*************************************************************************************

read.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/ERA5/concat/'
tmp.dir <- '/local_temp/ssobie/era5/'
if (!file.exists(tmp.dir))
   dir.create(tmp.dir,recursive=TRUE)

##-------------------------------------------------------------------------------------
##From Hourly Temperature
###'ssrd',
##var.list <- c('strd','tisr','fdir','uwind','vwind','dewpoint','tcc','lcc','rain')
##var.list <- 'wspd'
info.list <- list(pr=list(name='pr',type='pr',file='pr_hour_ERA5_BC_19800101-20181231.nc'),
                 sp=list(name='sp',type='sp',file='surface_pressure_hour_ERA5_BC_19800101-20181231.nc'),
                 ssrd=list(name='ssrd',type='ssrd',file='surface_solar_down_hour_ERA5_BC_19800101-20181231.nc'),
                 strd=list(name='strd',type='strd',file='surface_thermal_down_hour_ERA5_BC_19800101-20181231.nc'),
                 tisr=list(name='tisr',type='tisr',file='toa_insolation_hour_ERA5_BC_19800101-20181231.nc'),
                 fdir=list(name='fdir',type='fdir',file='total_sky_direct_hour_ERA5_BC_19800101-20181231.nc'),
                 uwind=list(name='uwind',type='uwind',file='uwind_hour_ERA5_BC_19800101-20181231.nc'),                 
                 vwind=list(name='vwind',type='vwind',file='vwind_hour_ERA5_BC_19800101-20181231.nc'),                 
                 wspd=list(name='wspd',type='wspd',file='wspd_hour_ERA5_BC_19800101-20181231.nc'),                 
                 dewpoint=list(name='dewpoint',type='dewpoint',file='dewpoint_hour_ERA5_BC_19800101-20181231.nc'),                 
                 tcc=list(name='tcc',type='tcc',file='total_cloud_cover_hour_ERA5_BC_19800101-20181231.nc'),                 
                 lcc=list(name='lcc',type='lcc',file='low_cloud_cover_hour_ERA5_BC_19800101-20181231.nc'),                 
                 rain=list(name='rain',type='rain',file='total_colum_rain_hour_ERA5_BC_19800101-20181231.nc'))

var.list <- 'pr'
var.fxn <- sum

for (var.name in var.list) {
  print(var.name)
  print(get_variable_units(var.name))
  var.info <- info.list[[var.name]]
  ##File for metadata attributes
  base.file <- var.info$file
  file.copy(from=paste0(read.dir,base.file),to=tmp.dir,overwrite=TRUE)

  calendar <- 'gregorian'
  hour.dates <- seq(from=as.PCICt('1980-01-01 00:00:00',cal=calendar),by='hour',to=as.PCICt('2018-12-31 23:00:00',cal=calendar))
  day.dates <- seq(from=as.PCICt('1980-01-01',cal=calendar),by='day',to=as.PCICt('2018-12-31',cal=calendar))
  years <- 1980:2018

  day.time <- era5_time_series('day',day.dates)

  daily.file <- make_daily_era5_netcdf(var.name,var.info$type,hour.dates=day.dates,base.file=base.file,
                                       dir=tmp.dir)

  dnc <- nc_open(paste0(tmp.dir,daily.file),write=TRUE)
  hnc <- nc_open(paste0(tmp.dir,base.file))
 
  for (year in years) {
     print(year)
     year.dates <- seq(from=as.PCICt(paste0(year,'-01-01 00:00:00'),cal=calendar),
                       by='hour',
                       to=as.PCICt(paste0(year,'-12-31 23:00:00'),cal=calendar))
     hix <- which(hour.dates %in% year.dates)
     hst <- head(hix,1)
     hen <- tail(hix,1)
     hct <- hen-hst+1
     year.data <- ncvar_get(hnc,var.name,start=c(1,1,hst),count=c(-1,-1,hct))
     daily.data <- make_daily_series(year.dates,year.data,var.fxn)

     yday.dates <- seq(from=as.PCICt(paste0(year,'-01-01 00:00:00'),cal=calendar),
                       by='day',
                       to=as.PCICt(paste0(year,'-12-31 23:00:00'),cal=calendar))
     yix <- which(day.dates %in% yday.dates)
     yst <- head(yix,1)
     yen <- tail(yix,1)
     yct <- yen-yst+1  

     ncvar_put(dnc,var.name,daily.data,start=c(1,1,yst),count=c(-1,-1,yct))
    
  }
  nc_close(dnc)
  nc_close(hnc)
  file.copy(from=paste0(tmp.dir,daily.file),to=read.dir,overwrite=TRUE)

}






