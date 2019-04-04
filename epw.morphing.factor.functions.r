##Script to generate shift and stretch factors for temperature
##from the BCCAQ2-TPS downscaled fields

library(ncdf4)
library(PCICt)
library(doParallel)
registerDoParallel(cores=4)
library(foreach)


source('/storage/home/ssobie/code/repos/crcm5/epw.belcher.functions.r',chdir=T)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##------------------------------------------------------------------------------
create.factor.file <- function(nc,var.name,var.units,
                               input,
                               tmp.dir,write.file) {

  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units

  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                     cal=time.calendar)
  time.values <- ncvar_get(nc,'time')
  time.series <- format(origin + time.values*86400,'%Y-%m-%d')
  years.ix <- grep('1995-',time.series)
  days <- time.values[years.ix]

  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  lon.atts <- ncatt_get(nc,'lon')
  lat.atts <- ncatt_get(nc,'lat')
  n.lon <- length(lon)
  n.lat <- length(lat)

  global.atts <- ncatt_get(nc,varid=0)
  ##--------------------------------------------------------------
  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units, days,
                      unlim=TRUE, calendar=time.calendar)
  var.geog <- ncvar_def(var.name, units=var.units, dim=list(x.geog, y.geog, t.geog),
                        missval=-32768)
  file.nc <- nc_create(paste(tmp.dir,write.file,sep=''), var.geog)
  ##File Attributes
  global.names <- names(global.atts)
  for (g in 1:length(global.atts))
    ncatt_put(file.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])

  ##Time attributes
  ncatt_put(file.nc,varid='time',attname='units',attval=time.units)
  ncatt_put(file.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='calendar',attval=time.calendar)

  lon.names <- names(lon.atts)
  for (j in 1:length(lon.atts))
    ncatt_put(file.nc,varid='lon',attname=lon.names[j],attval=lon.atts[[j]])

  lat.names <- names(lat.atts)
  for (j in 1:length(lat.atts))
    ncatt_put(file.nc,varid='lat',attname=lat.names[j],attval=lat.atts[[j]])

  ##Climdex Attributes
  ncatt_put(file.nc,varid=var.name,attname='units',attval=var.units)
  ncatt_put(file.nc,varid=var.name,attname='_FillValue',attval=-32768)
  ncatt_put(file.nc,varid=var.name,attname='standard_name',attval=var.name)
  ncatt_put(file.nc,varid=var.name,attname='long_name',attval=var.name)

  ncatt_put(file.nc,varid=0,attname='history',attval='')
  ncvar_put(file.nc,var.name,input)
  ncvar_put(file.nc,'lon',lon)
  ncvar_put(file.nc,'lat',lat)

  nc_close(file.nc)

}

##------------------------------------------------------------------------------
make_average_series <- function(series,fac,method,rlen=1,agg.fxn) {

   agg.daily <- tapply(series,fac,agg.fxn)
   if (method=='roll') {
     agg.daily <- rollmean(agg.data,as.numeric(rlen),fill='extend')
   }
  return(agg.daily)
}

##------------------------------------------------------------------------------


fix.time.series <- function(nc,gcm,interval,method) {
   time.atts <- ncatt_get(nc,'time')
   time.calendar <- time.atts$calendar
   time.units <- time.atts$units
   time.values <- ncvar_get(nc,'time')
   origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                            cal=time.calendar)
   time.series <- origin.pcict + time.values*86400
   years <- format(time.series,'%Y')
   yrs <- strsplit(interval,'-')[[1]]

   feb.flag <- grep('-02-29',time.series)

   new.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                          cal='365_day')
                           
   new.time0 <- (as.PCICt(format(time.series[1],'%Y-%m-%d'),cal='365_day') - new.origin)/86400
   if (length(feb.flag)==0) {
      new.values <- seq(as.numeric(new.time0),by=1,length.out=length(time.values))
   } else {
      new.values <- seq(as.numeric(new.time0),by=1,length.out=length(time.values[-feb.flag]))
   }

   new.series <- new.origin + new.values*86400

   years <- format(new.series,'%Y')
   st <- head(grep(yrs[1],years),1)
   en <- tail(grep(yrs[2],years),1)
   if (grepl('HadGEM',gcm) & yrs[2]=='2100') {
     en <- length(years)
   }
   cnt <- en-st+1

   factor <- switch(method,
                    daily='%m-%d',
                    monthly='%m',
                    roll='%m-%d',
                    seasonal='%m')
   fac <- as.factor(format(new.series[st:en],factor))

   rv <- list(time=new.series[st:en],fac=fac,st=st,en=en,cnt=cnt)
   return(rv)
}




##----------------------------------------------------------------------------------------------
separate.into.list <- function(nc,var.name,j,time,time.bnds) {

    data.raw <- ncvar_get(nc,var.name,start=c(1,j,time.bnds$st),count=c(-1,1,time.bnds$cnt))
    feb.flag <- grep('-02-29',time.bnds$time)
    if (length(feb.flag!=0)) {
      data.subset <- data.raw[-feb.flag]
    } else {
      data.subset <- data.raw
    }
    data.list <- lapply(seq_len(nrow(data.subset)), function(k) data.subset[k,])
    rm(data.subset)
    rm(data.raw)
    return(data.list)
}

##----------------------------------------------------------------------------------------------
##Calculate morphing factors

daily.aggregate <- function(nc,var.name,
                            gcm,interval,
                            method,rlen,agg.fxn) {

   time <- netcdf.calendar(nc)
   n.lat <- nc$dim$lat$len ##Latitude Length
   n.lon <- nc$dim$lon$len ##Longitude Length
   n.time <- length(time)

   time.bnds <- fix.time.series(nc,gcm,interval,method)
   sub.time <- time.bnds$time
   sub.fac <- time.bnds$fac

   agg.array <- array(NA,c(n.lon,n.lat,length(levels(sub.fac))))

   ##Iterate along latitude indices
   for (j in 1:n.lat) {
      print(paste0('Latitude: ',j,' of ',n.lat))

      data.list <- separate.into.list(nc,var.name,j,time,time.bnds)

      na.flag <- unlist(lapply(data.list,function(x){any(is.na(x))}))
      data.sub.list <- data.list[!na.flag]

      data.test <- make_average_series(series=data.sub.list[[1]],
                                         fac=sub.fac,method=method,rlen=rlen,agg.fxn=agg.fxn)           

      data.result <- rep(list(rep(NA,length(levels(sub.fac)))),n.lon)

      ##Compute the aggregate values (default to daily)
      data.agg <- foreach(
                      data=data.sub.list,
                      .export=c('make_average_series','sub.fac','method','rlen','agg.fxn')
                    ) %dopar% {
                      objects <- make_average_series(series=data,
                                      fac=sub.fac,method=method,rlen=rlen,agg.fxn=agg.fxn)
                              }
      data.result[!na.flag] <- data.agg
      ncol <- length(data.result[[1]])
      data.matrix <- matrix(unlist(data.result),nrow=n.lon,ncol=ncol,byrow=TRUE)
      agg.array[,j,] <- data.matrix
   }
   return(agg.array)
}

##----------------------------------------------------------------------------------------------

