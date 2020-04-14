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
                    
  rv <- hourly.time

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
                tasskew='',sphour='Pa',uwind='m s-1',vwind='m s-1',dewpoint='degC',
                tcc='fraction',lcc='fraction',rain='kg m-2',
                ssrd='kJ m-2',strd='kJ m-2',tisr='kJ m-2',fdir='kJ m-2')             
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
  sphour.atts <- list(long_name = "Surface pressure",
                       cell_methods = "time: skewness",units=get_variable_units('tasskew'))
  uwind.atts <- list(long_name = "10 metre U wind component",
                     units=get_variable_units('uwind'))
  vwind.atts <- list(long_name = "10 metre V wind component",
                     units=get_variable_units('vwind'))
  dewpoint.atts <- list(long_name = "2 metre dewpoint temperature",
                        units=get_variable_units('dewpoint'))
  totalcloud.atts <- list(long_name = "Total cloud cover",
                        units=get_variable_units('totalcloud'))
  lowcloud.atts <- list(long_name = "Low cloud cover",
                        units=get_variable_units('lowcloud'))
  totalrain.atts <- list(long_name = "Total Column Rain Water",
                        units=get_variable_units('totalrain'))
  ssrd.atts <- list(long_name = "Surface solar radiation downwards",
                        units=get_variable_units('ssrd'))
  strd.atts <- list(long_name = "Surface thermal radiation downwards",
                        units=get_variable_units('strd'))
  tisr.atts <- list(long_name = "TOA incident solar radiation",
                        units=get_variable_units('tisr'))
  fdir.atts <- list(long_name = "Total sky direct solar radiation at surface",
                        units=get_variable_units('fdir'))

  var.atts <- switch(var.name,
                     prhour=pr.hour.atts,
                     prday=pr.day.atts,
                     tasmax=tasmax.atts,
                     tasmin=tasmin.atts,
                     tasday=tasday.atts,
                     tashour=tashour.atts,
                     tasrange=tasrange.atts,
                     tasskew=tasskew.atts,
                     uwind=uwind.atts,
                     vwind=vwind.atts,
                     dewpoint=dewpoint.atts,
                     tcc=totalcloud.atts,
                     lcc=lowcloud.atts,
                     rain=totalrain.atts,
                     ssrd=ssrd.atts,
                     strd=strd.atts,
                     tisr=tisr.atts,
                     fdir=fdir.atts)
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
  dewpoint.atts <- list(standard_name = "2 metre dewpoint temperature",
                   missing_value = 1.e+20)
  totalcloud.atts <- list(standard_name = "Cloud cover fraction",
                   missing_value = 1.e+20)
  lowcloud.atts <- list(standard_name = "Cloud cover fraction",
                   missing_value = 1.e+20)
  totalrain.atts <- list(standard_name = "Total Column Rain Water",
                   missing_value = 1.e+20)
  ssrd.atts <- list(standard_name = "Surface solar radiation downwards",
                    missing_value = 1.e+20)                        
  strd.atts <- list(standard_name = "Surface thermal radiation downwards",
                   missing_value = 1.e+20)
  tisr.atts <- list(standard_name = "TOA incident solar radiation",
                   missing_value = 1.e+20)
  fdir.atts <- list(standard_name = "Total sky direct solar radiation at surface",
                   missing_value = 1.e+20)


  var.atts <- switch(var.name,
                     pr=pr.atts,
                     tas=tas.atts,
                     sp=sp.atts,
                     uwind=uwind.atts,
                     vwind=vwind.atts,
                     dewpoint=dewpoint.atts,
                     totalcloud=totalcloud.atts,
                     lowcloud=lowcloud.atts,
                     totalrain=totalrain.atts,
                     ssrd=ssrd.atts,
                     strd=strd.atts,
                     tisr=tisr.atts,
                     fdir=fdir.atts)

  rv <- list(lon=lon.atts,
             lat=lat.atts,
             var=var.atts)
  return(rv)
}

add_attributes_ncdf <- function(var.info, time, nc) {

  standard.atts <- get_standard_atts(var.info$name)
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
    ncatt_put(nc,varid=var.info$name,attname=var.names[j],attval=standard.atts$var[[j]])
  print('Variable names')
  var.names <- names(variable.atts$var)
  for (j in 1:length(variable.atts$var))
    ncatt_put(nc,varid=var.info$name,attname=var.names[j],attval=variable.atts$var[[j]])
 
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

  var.geog <- ncvar_def(var.info$name, units=get_variable_units(var.info$name),
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
concatenate_series <- function(var.name,month,vnc,var.info,input.data) {
   print(var.name)
   all.dates <- format(var.info$dates,'%Y-%m')
   fmon <- gsub('_','-',month)
   yix <- which(all.dates %in% fmon)
   yst <- head(yix,1)
   yen <- tail(yix,1)
   yct <- yen-yst+1
   ncvar_put(vnc,var.info$name,round(input.data/1000,3),start=c(1,1,yst),count=c(-1,-1,yct))
}
 
##-------------------------------------------------------------------------------------


##*************************************************************************************

read.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/ERA5/'
tmp.dir <- '/local_temp/ssobie/era55/'
if (!file.exists(tmp.dir))
   dir.create(tmp.dir,recursive=TRUE)

##-------------------------------------------------------------------------------------
##From Hourly Temperature

##File for metadata attributes
base.file <- 'radiation_variables_hour_ERA5_BC_1995_01.nc'
file.copy(from=paste0(read.dir,base.file),to=tmp.dir,overwrite=TRUE)

calendar <- 'gregorian'
hour.dates <- seq(from=as.PCICt('1980-01-01 00:00:00',cal=calendar),by='hour',to=as.PCICt('2018-12-31 23:00:00',cal=calendar))

hour.time <- era5_time_series('hour',hour.dates)

var.list <- c('ssrd','strd','tisr','fdir')

info.list <- list(ssrd=list(name='ssrd',filevar='ssrd',type='rad',freq='day',dates=hour.dates,
                                infile='radiation_variables_hour_ERA5_BC', file='surface_solar_down_day_ERA5_BC_19800101-20181231.nc'),
                  strd=list(name='strd',filevar='strd',type='rad',freq='day',dates=hour.dates,
                                infile='radiation_variables_hour_ERA5_BC',file='surface_thermal_down_day_ERA5_BC_19800101-20181231.nc'),
                  tisr=list(name='tisr',filevar='tisr',type='rad',freq='day',dates=hour.dates,
                                infile='radiation_variables_hour_ERA5_BC',file='toa_insolation_day_ERA5_BC_19800101-20181231.nc'),
                  fdir=list(name='fdir',filevar='fdir',type='rad',freq='day',dates=hour.dates,
                                infile='radiation_variables_hour_ERA5_BC',file='total_sky_direct_day_ERA5_BC_19800101-20181231.nc'))

##----------------------------
##Create derived files
ncs.list <- vector(mode='list',length=length(var.list))
names(ncs.list) <- var.list
for (var.name in var.list) {
   print(paste0('Making file for ',var.name))
   var.info <- info.list[[var.name]]
   make_era5_netcdf(var.info,base.file=base.file,
                    dir=tmp.dir)
   ncs.list[[var.name]] <- nc_open(paste0(tmp.dir,var.info$file),write=TRUE)
 
}

##---------------------------

  for (month in months) {
    print(month)
    month.file <- paste0('radiation_variables_hour_ERA5_BC_',month,'.nc')
    file.copy(from=paste0(read.dir,month.file),to=tmp.dir,overwrite=TRUE)

    ync <- nc_open(paste0(tmp.dir,month.file))
    ##Invert tas to fix latitude order
    n.lat <- ync$dim$latitude$len

    ##Concatenate Hourly temperatures
    ssrd.raw <- ncvar_get(ync,'ssrd')
    ssrd.data <- ssrd.raw[,(n.lat:1),]
    ssrd.info <- info.list[['ssrd']]
    concatenate_series('ssrd',month,ncs.list[['ssrd']],ssrd.info,ssrd.data)

    strd.raw <- ncvar_get(ync,'strd')
    strd.data <- strd.raw[,(n.lat:1),]
    strd.info <- info.list[['strd']]
    concatenate_series('strd',month,ncs.list[['strd']],strd.info,strd.data)

    tisr.raw <- ncvar_get(ync,'tisr')
    tisr.data <- tisr.raw[,(n.lat:1),]
    tisr.info <- info.list[['tisr']]
    concatenate_series('tisr',month,ncs.list[['tisr']],tisr.info,tisr.data)

    fdir.raw <- ncvar_get(ync,'fdir')
    fdir.data <- fdir.raw[,(n.lat:1),]
    fdir.info <- info.list[['fdir']]
    concatenate_series('fdir',month,ncs.list[['fdir']],fdir.info,fdir.data)

    nc_close(ync)
    file.remove(paste0(tmp.dir,month.file))
  }


for (var.name in var.list) {
  nc_close(ncs.list[[var.name]])
  print(info.list[[var.name]][['file']])
  file.copy(from=paste0(tmp.dir,info.list[[var.name]][['file']]),to=paste0(read.dir,'concat/'),overwrite=TRUE)
}




