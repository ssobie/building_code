##Script to calculate dewpoint temperatures

library(ncdf4)
library(PCICt)
library(doParallel)
registerDoParallel(cores=4)
library(foreach)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
##source('/storage/home/ssobie/code/repos/crcm5/add.crcm5.metadata.r',chdir=T)

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


##Wet Bulb Temperature
temp.wet.bulb <- function(tas,temp.dew,pas,sp.hum) {

     ##Mixing Ratio
     mix.ratio <- sp.hum / (1 - sp.hum)

     ##Specific Heat Capacity
     heat.cap <- cpd * (1 + 0.84*mix.ratio)

     ##Saturation Vapour Pressure
     vp.sat <- sat.vape * exp( (lh.vape/R.vapor) * (1/T.zero - 1/tas))

     ##Saturation Mixing Ratio
     sat.mix.ratio <- (epsilon * vp.sat) / (pas - vp.sat)

     ##Saturated Adiabatic Lapse Rate
     a.lapse.sat <- (grav/heat.cap) * (1 + (sat.mix.ratio * lh.vape)/(R.dry * tas)) /
                                      (1 + (lh.vape^2 * sat.mix.ratio * epsilon)/(heat.cap * R.dry * tas^2))

     ##Lifting Condensation Level Height
     lift.height <- a.lapse * (tas - temp.dew)

     ##Air temperature at the lifting condensation level
     temp.lcl <- tas - (d.lapse * lift.height)

     ##Wet Bulb Temperature
     temp.WB <- temp.lcl + (a.lapse.sat * lift.height)
     return(temp.WB)
}



##----------------------------------------------------------------------------------------------

make.new.wbt.file <- function(gcm,rcp,run,
                              tas.file,
                              tmp.base) {
  new.var <- 'wetbulb'
  nc <- nc_open(paste0(tmp.base,tas.file),write=FALSE)
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*86400

  yst <-  gsub('-','',format(head(time.series,1),'%Y-%m-%d'))
  yen <-  gsub('-','',format(tail(time.series,1),'%Y-%m-%d'))

  write.file <- paste0('wetbulb_day_',gcm,'_historical+',scenario,'_',run,'_19500101-21001231.nc')
  ##Attributes to retain
  lon <- ncvar_get(nc,'lon')
  lon_bnds <- ncvar_get(nc,'lon_bnds')
  lat <- ncvar_get(nc,'lat')
  lat_bnds <- ncvar_get(nc,'lat_bnds')
  var.atts <- ncatt_get(nc,'tas')

  time_bnds <- ncvar_get(nc,'time_bnds')
  n.lon <- length(lon)
  n.lat <- length(lat)

  ##--------------------------------------------------------------
  ##Create new netcdf file

  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units, time.values,
                      unlim=FALSE, calendar=time.calendar)
  var.geog <- ncvar_def(new.var, units=var.atts$units, dim=list(x.geog, y.geog,t.geog),
                        missval=1.e+20)

  hist.nc <- nc_create(paste(tmp.base,write.file,sep=''),var.geog)

  ##Add Attributes
  print('Time atts')
  time.atts <- ncatt_get(nc,'time')
  time.names <- names(time.atts)
  for (t in 1:length(time.atts))
    ncatt_put(hist.nc,varid='time',attname=time.names[t],attval=time.atts[[t]])

  lon.atts <- ncatt_get(nc,'lon')
  print(lon.atts)
  lon.names <- names(lon.atts)
  print('Lon names')
  for (j in 1:4) ##length(lon.atts))
    ncatt_put(hist.nc,varid='lon',attname=lon.names[j],attval=lon.atts[[j]])

  lat.atts <- ncatt_get(nc,'lat')
  print('Lat names')
  lat.names <- names(lat.atts)
  for (j in 1:4) ##length(lat.atts))
    ncatt_put(hist.nc,varid='lat',attname=lat.names[j],attval=lat.atts[[j]])

  print('Global Atts')
  global.atts <- ncatt_get(nc,0)
  global.names <- names(global.atts)
  for (g in 1:length(global.atts))
    ncatt_put(hist.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])

  print('Variable Atts')
  dewpoint.atts <- list(standard_name = "wetbulb_temperature",
                   long_name = "Near-Surface Wetbulb Temperature",
                   units = "degC",
                   missing_value =1.e+20,
                   cell_methods = "time: mean")
                   
  varnames <- names(dewpoint.atts)
  for (j in 1:length(dewpoint.atts))
    ncatt_put(hist.nc,varid=new.var,attname=varnames[j],attval=dewpoint.atts[[j]])

  ##Clear extraneous history
  ncatt_put(hist.nc,varid=0,attname='history',attval='')

  nc_close(hist.nc)
  nc_close(nc)

  return(write.file)

}



##-----------------------------------------------------------------------------------------

calculate.wetbulb.temp <- function(tas.file,dwpt.file,psl.file,huss.file,wetbulb.file,tmp.base,offset) {

  tas.nc <- nc_open(paste0(tmp.base,tas.file))
  dwpt.nc <- nc_open(paste0(tmp.base,dwpt.file))
  psl.nc <- nc_open(paste0(tmp.base,psl.file))
  huss.nc <- nc_open(paste0(tmp.base,huss.file))
  wt.nc <- nc_open(paste0(tmp.base,wetbulb.file),write=TRUE)

  time <- netcdf.calendar(tas.nc)
  yrs <- levels(as.factor(format(time,'%Y')))
  n.lat <- tas.nc$dim$lat$len ##Latitude Length
  n.lon <- tas.nc$dim$lon$len ##Longitude Length
  n.time <- tas.nc$dim$time$len ##Longitude Length

  for (j in 1:n.lat) { 
    print(paste0('Latitude: ',j,' of ',n.lat))
    tas.subset <- ncvar_get(tas.nc,'tas',start=c(1,j,1),count=c(-1,1,-1))
    tas.list <- lapply(seq_len(nrow(tas.subset)), function(k) tas.subset[k,])
    rm(tas.subset)

    dwpt.subset <- ncvar_get(dwpt.nc,'dewpoint',start=c(1,j+offset,1),count=c(-1,1,-1))
    dwpt.list <- lapply(seq_len(nrow(dwpt.subset)), function(k) dwpt.subset[k,])
    rm(dwpt.subset)

    psl.subset <- ncvar_get(psl.nc,'psl',start=c(1,j+offset,1),count=c(-1,1,-1))
    psl.list <- lapply(seq_len(nrow(psl.subset)), function(k) psl.subset[k,])
    rm(psl.subset)

    huss.subset <- ncvar_get(huss.nc,'huss',start=c(1,j+offset,1),count=c(-1,1,-1))
    huss.list <- lapply(seq_len(nrow(huss.subset)), function(k) huss.subset[k,])
    rm(huss.subset)

    wetbulbs <- foreach(
                       tas=tas.list,                 
                       dwpt=dwpt.list,                 
                       psl=psl.list,                 
                       huss=huss.list,                 
                       .export=c('temp.wet.bulb')
                             ) %dopar% {
                                objects <- temp.wet.bulb(tas=tas,temp.dew=dwpt,pas=psl,sp.hum=huss)
                                    }
    wt.matrix <- matrix(unlist(wetbulbs),nrow=n.lon,ncol=n.time,byrow=T)
    ncvar_put(wt.nc,varid='wetbulb',vals=wt.matrix,
                      start=c(1,j,1),count=c(-1,1,n.time))
  }
  nc_close(tas.nc)
  nc_close(dwpt.nc)
  nc_close(psl.nc)
  nc_close(huss.nc)
}




##----------------------------------------------------------------------------------------------
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }


tmp.dir <- tmpdir
print(tmp.dir)
print(tmpdir)
if (!file.exists(tmpdir)) {
  dir.create(tmpdir,recursive=T)
}

run <- 'r3i1p1'

read.dir <- paste0('/storage/data/climate/downscale/CMIP5/building_code/',gcm,'/')
dwpt.file <- paste0('dewpoint_day_',gcm,'_historical+',scenario,'_',run,'_19500101-21001231.nc')
psl.file <- paste0('psl_day_',gcm,'_historical+',scenario,'_',run,'_19500101-21001231.nc')
huss.file <- paste0('huss_day_',gcm,'_historical+',scenario,'_',run,'_19500101-21001231.nc')

tas.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/CMIP5/',gcm,'/')
tas.file <- paste0('tas_day_',gcm,'_historical+rcp85_r1i1p1_19500101-21001231.nc')

  nc.tas <- nc_open(paste0(tas.dir,tas.file))
  tas.lat <- ncvar_get(nc.tas,'lat')

  nc.psl <- nc_open(paste0(read.dir,psl.file))
  psl.lat <- ncvar_get(nc.psl,'lat')
  offset <- which(psl.lat %in% tas.lat[1]) - 1
  if (offset < 1) 
    offset <- 0

print(offset)

nc_close(nc.psl)
nc_close(nc.tas)

if (1==1) {
print('Copying TAS file')
file.copy(from=paste0(tas.dir,tas.file),to=tmp.dir,overwrite=T)
print('Copying DWPT file')
file.copy(from=paste0(read.dir,dwpt.file),to=tmp.dir,overwrite=T)
print('Copying PSL file')
file.copy(from=paste0(read.dir,psl.file),to=tmp.dir,overwrite=T)
print('Copying HUSS file')
file.copy(from=paste0(read.dir,huss.file),to=tmp.dir,overwrite=T)

wetbulb.file <- make.new.wbt.file(gcm,scenario,run,
                                  tas.file,
                                  tmp.dir)
calculate.wetbulb.temp(tas.file,dwpt.file,psl.file,huss.file,wetbulb.file,tmp.dir,offset)

file.copy(from=paste0(tmp.dir,wetbulb.file),to=read.dir,overwrite=T)
print(wetbulb.file)
file.remove(paste0(tmp.dir,tas.file))
print('Done with TAS')
file.remove(paste0(tmp.dir,dwpt.file))
print('Done with Dewpoint')
file.remove(paste0(tmp.dir,psl.file))
print('Done with PSL')
file.remove(paste0(tmp.dir,huss.file))
print('Done with HUSS')
file.remove(paste0(tmp.dir,wetbulb.file))
print('Done with Wetbulb')

}


