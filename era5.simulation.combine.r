##Script to reorder the ERA5 files to facilitate a 
##sliding window calibration option for ClimDown


library(ncdf4)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R')

##-------------------------------------------------------------------------------------
##Calibration intervals:
##1980-2010 -> Sim 2011-2018
##1990-2018 -> Sim 1980-1989
##1980-1989 and 2000-2018 -> Sim 1990-1999
##1980-1999 and 2010-2018 -> Sim 2000-2009

##-------------------------------------------------------------------------------------
##Uses Gregorian Calendar
era5_time_series <- function(cstart1,cend1,cstart2=NULL,cend2=NULL) {

  day.dates <- seq(from=as.Date('1980-01-01'),by='day',to=as.Date('2018-12-31'))
  cst1 <- head(grep(cstart1,day.dates),1)
  cen1 <- tail(grep(cend1,day.dates),1)
  sub1 <- day.dates[cst1:cen1]
  sub.dates <- sub1
  if (!is.null(cstart2)) {
     cst2 <- head(grep(cstart2,day.dates),1)
     cen2 <- tail(grep(cend2,day.dates),1)
     sub2 <- day.dates[cst2:cen2]
     sub.dates <- sort(c(sub1,sub2))
  } 

  ix <- day.dates %in% sub.dates
  print(range(day.dates[ix]))
  print(range(day.dates[!ix]))
  return(ix)
}


read.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/ERA5/concat/'
write.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/ERA5/concat/sim/'
tmp.dir <- '/local_temp/ssobie/era5/'
if (!file.exists(tmp.dir))
   dir.create(tmp.dir,recursive=TRUE)

var.list <- c('tas','tasmax','tasmin','tasrange','tasskew')
cal.list <- list(list(cstart1='1980',cend1='2009',cstart2=NULL,cend2=NULL,sim='20110101-20181231'),
                 list(cstart1='1990',cend1='2018',cstart2=NULL,cend2=NULL,sim='19800101-19891231'),
                 list(cstart1='1980',cend1='1989',cstart2='2000',cend2='2018',sim='19900101-19991231'),
                 list(cstart1='1980',cend1='1999',cstart2='2010',cend2='2018',sim='20000101-20091231'))

for (var.name in var.list) {
    print(var.name)
    base.file <- paste0(var.name,'_day_ERA5_BC_19800101-20181231.nc')   
    write.file <- paste0(var.name,'_QDM_day_ERA5_BC_19800101-20181231.nc')   
    file.copy(from=paste0(read.dir,base.file),to=paste0(tmp.dir,write.file),overwrite=TRUE)
    nc <- nc_open(paste0(tmp.dir,write.file),write=TRUE)
    data <- ncvar_get(nc,var.name)
    time <- netcdf.calendar(nc)
    sim.data <- data*0
    for (cal in cal.list) {
       print(cal$sim)

       sim.file  <- paste0(var.name,'_day_QDM_ERA5_',cal$sim,'.nc')   
       file.copy(from=paste0(read.dir,'sim/',sim.file),to=tmp.dir,overwrite=TRUE)

       snc <- nc_open(paste0(tmp.dir,sim.file))

       cal.ix <- era5_time_series(cstart1=cal$cstart1,cend1=cal$cend1,cstart2=cal$cstart2,cend2=cal$cend2)

       cal.len <- sum(cal.ix)
       all.len <- length(cal.ix)

       sim.data[,,!cal.ix] <- data[,,(cal.len+1):all.len]
       nc_close(snc)
       file.remove(paste0(tmp.dir,sim.file))
    }
    ncvar_put(nc,var.name,sim.data)
    nc_close(nc)
    file.copy(from=paste0(tmp.dir,write.file),to=write.dir,overwrite=TRUE)
}    