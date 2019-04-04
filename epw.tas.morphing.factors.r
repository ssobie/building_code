##Script to generate shift and stretch factors for temperature
##from the BCCAQ2-TPS downscaled fields

library(ncdf4)
library(PCICt)
library(doParallel)
registerDoParallel(cores=4)
library(foreach)

source('/storage/home/ssobie/code/repos/building_code/epw.morphing.factor.functions.r')
##source('/storage/home/ssobie/code/repos/crcm5/epw.belcher.functions.r',chdir=T)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##**************************************************************************************
##**************************************************************************************

scenario <- 'TPS'
method <- 'daily'
rlen <- '21'
agg.fxn <- mean
##gcm <- 'ACCESS1-0'

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

if (method!='roll') { 
  rlen <- ''
}

tas.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/'
epw.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/epw_factors/'

tmp.dir <- '/local_temp/ssobie/epw/'
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}

##--------------------------------------------------------------------------------------------
##Copy temperature files to temporary directory

##dir.create(paste0(tmp.dir,'/',gcm,'/'))
scen.files <- list.files(path=paste0(tas.dir,gcm,'/'),pattern=scenario)
tasmax.file <- scen.files[grep('tasmax_BCCAQ2',scen.files)]
file.copy(from=paste0(tas.dir,gcm,'/',tasmax.file),to=paste0(tmp.dir,'/'),overwrite=TRUE)
tasmin.file <- scen.files[grep('tasmin_BCCAQ2',scen.files)]
file.copy(from=paste0(tas.dir,gcm,'/',tasmin.file),to=paste0(tmp.dir,'/'),overwrite=TRUE)

##------------------------------------------------------------------------------


tasmax.nc <- nc_open(paste0(tmp.dir,'/',tasmax.file))
tasmin.nc <- nc_open(paste0(tmp.dir,'/',tasmin.file))

##
tasmax.1980s <- daily.aggregate(tasmax.nc,'tasmax',
                        gcm=gcm,interval='1971-2000',
                        method=method,rlen=rlen,agg.fxn=agg.fxn)

tasmin.1980s <- daily.aggregate(tasmin.nc,'tasmin',
                        gcm=gcm,interval='1971-2000',
                        method=method,rlen=rlen,agg.fxn=agg.fxn)
tas.1980s <- (tasmax.1980s + tasmin.1980s)/2

intervals <- c('2011-2040','2041-2070','2071-2100')

for (interval in intervals) {

  tasmax.proj <- daily.aggregate(tasmax.nc,'tasmax',
                        gcm=gcm,interval=interval,
                        method=method,rlen=rlen,agg.fxn=agg.fxn)

  tasmin.proj <- daily.aggregate(tasmin.nc,'tasmin',
                        gcm=gcm,interval=interval,
                        method=method,rlen=rlen,agg.fxn=agg.fxn)
  tas.proj <- (tasmax.proj + tasmin.proj)/2

  delta.tasmax <- tasmax.proj - tasmax.1980s

  deltx.file <- paste0('delta_tasmax_',gcm,'_1971-2000_',interval,'.nc')
  create.factor.file(tasmax.nc,'tasmax','degC',
                     delta.tasmax,
                     tmp.dir,deltx.file)
  file.copy(from=paste0(tmp.dir,deltx.file),to=epw.dir,overwrite=TRUE)

  delta.tasmin <- tasmin.proj - tasmin.1980s
  deltn.file <- paste0('delta_tasmin_',gcm,'_1971-2000_',interval,'.nc')
  create.factor.file(tasmin.nc,'tasmin','degC',
                   delta.tasmin,
                   tmp.dir,deltn.file)
  file.copy(from=paste0(tmp.dir,deltn.file),to=epw.dir,overwrite=TRUE)

  delta.tas <- tas.proj - tas.1980s
  delts.file <- paste0('delta_tas_',gcm,'_1971-2000_',proj,'.nc')
  create.factor.file(tasmin.nc,'tas','degC',
                     delta.tas,
                     tmp.dir,delts.file)
  file.copy(from=paste0(tmp.dir,delts.file),to=epw.dir,overwrite=TRUE)

  alpha.numer <- (delta.tasmax - delta.tasmin)
  alpha.file <- paste0('alpha_tasmax_tasmin_',gcm,'_1971-2000_',interval,'.nc')
  create.factor.file(tasmax.nc,'alpha_tas','decC',
                   alpha.numer,
                   tmp.dir,alpha.file)
  file.copy(from=paste0(tmp.dir,alpha.file),to=epw.dir,overwrite=TRUE)

  file.remove(paste0(tmp.dir,alpha.file))
  file.remove(paste0(tmp.dir,deltx.file))
  file.remove(paste0(tmp.dir,deltn.file))
  file.remove(paste0(tmp.dir,delts.file))

}

nc_close(tasmax.nc)
nc_close(tasmin.nc)


file.remove(paste0(tmp.dir,tasmax.file))
file.remove(paste0(tmp.dir,tasmin.file))
