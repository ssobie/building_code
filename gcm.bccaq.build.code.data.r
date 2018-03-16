##Contains functions to compute Building Code values from the BC
##Building Code and the American ASHRAE Handbook

library(ncdf4)
library(PCICt)

get.time.series <- function(nc) {

  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units

  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  values <- ncvar_get(nc,'time')
  series <- format(origin + (values)*86400,'%Y-%m-%d')
  return(series)
}

read.data.subset <- function(data,dates,interval) {

  bnds <- strsplit(interval,'-')[[1]]
  st <- head(grep(bnds[1],dates),1)
  en <- tail(grep(bnds[2],dates),1)
  if (length(en)==0) {
    en <- length(dates)
  }
  data.subset <- data[st:en]
  dates.subset <- dates[st:en]
  rv <- data.subset
  return(rv)
}

get.gcm.data <- function(var.name,gcm,gcm.file,lon.c,lat.c) {

  print(gcm.file)
  nc <- nc_open(gcm.file)
  lon <- ncvar_get(nc,'lon')
  lon <- ((lon + 180) %% 360) - 180
  lat <- ncvar_get(nc,'lat')
  
  lon.ix <- which.min(abs(lon.c-lon))
  lat.ix <- which.min(abs(lat.c-lat))
  gcm.data <- ncvar_get(nc,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  gcm.time <- get.time.series(nc)
  gcm.subset <- read.data.subset(gcm.data,gcm.time,'1951-2100')
  time.subset <- read.data.subset(gcm.time,gcm.time,'1951-2100')
  nc_close(nc)
  if (grepl(var.name,'(tas|tasmax|tasmin)'))
     gcm.subset <- gcm.subset - 273
  if (grepl(var.name,'(pr)')) {
     gcm.subset <- gcm.subset * 86400
     gcm.subset[gcm.subset <0] <- 0
  }
  if (grepl(var.name,'(psl)'))
     gcm.subset <- gcm.subset/1000

  rv <- list(data=gcm.subset,time=time.subset)

  return(rv)
}

get.snow.data <- function(var.name,gcm,scenario,data.dir,lon.c,lat.c) {

   gcm.files <- list.files(path=data.dir,
                           pattern=paste(var.name,'_BCCAQ-PRISM_',gcm,'_',scenario,sep=''),full.name=TRUE)

   print(gcm.files)
   past.file <- gcm.files[grep('1951-2100',gcm.files)]
   nc.past <- nc_open(past.file)

   lon <- ncvar_get(nc.past,'lon')
   lon <- ((lon + 180) %% 360) - 180
   lat <- ncvar_get(nc.past,'lat')
  
   lon.ix <- which.min(abs(lon.c-lon))
   lat.ix <- which.min(abs(lat.c-lat))

   past.subset <- ncvar_get(nc.past,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
   past.time <- get.time.series(nc.past)
   nc_close(nc.past)


   rv <- list(data=past.subset,time=past.time)
   return(rv)
}


get.bccaq.data <- function(var.name,gcm,scenario,data.dir,lon.c,lat.c) {

   gcm.files <- list.files(path=data.dir,
                           pattern=paste(var.name,'_day_BCCAQ_',gcm,'_',scenario,sep=''),full.name=TRUE)

   print(gcm.files)
   past.file <- gcm.files[grep('1951-2000',gcm.files)]
   proj.file <- gcm.files[grep('2001-2100',gcm.files)]
   nc.past <- nc_open(past.file)
   nc.proj <- nc_open(proj.file)

   lon <- ncvar_get(nc.past,'lon')
   lon <- ((lon + 180) %% 360) - 180
   lat <- ncvar_get(nc.past,'lat')
  
   lon.ix <- which.min(abs(lon.c-lon))
   lat.ix <- which.min(abs(lat.c-lat))


   past.subset <- ncvar_get(nc.past,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
   past.time <- get.time.series(nc.past)
   proj.subset <- ncvar_get(nc.proj,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
   proj.time <- get.time.series(nc.proj)
   comb.subset <- c(past.subset,proj.subset)
   time.subset <- as.Date(c(past.time,proj.time))

   nc_close(nc.past)
   nc_close(nc.proj)

   rv <- list(data=comb.subset,time=time.subset)
   return(rv)
}


gather.gcm.data <- function(gcm,lon.c,lat.c,scenario) {
  gcm.dir <- paste('/storage/data/climate/downscale/CMIP5/building_code/',gcm,sep='')

  rhs.file <- list.files(path=gcm.dir,pattern=paste0('rhs_day_',gcm),full.name=TRUE)
  rhs.data <- get.gcm.data(var.name='rhs',gcm=gcm,gcm.file=rhs.file,lon.c,lat.c)

  rsds.file <- list.files(path=gcm.dir,pattern=paste0('rsds_day_',gcm),full.name=TRUE)
  rsds.data <- get.gcm.data(var.name='rsds',gcm=gcm,gcm.file=rsds.file,lon.c,lat.c)

  clt.file <- list.files(path=gcm.dir,pattern=paste0('clt_day_',gcm),full.name=TRUE)
  clt.data <- get.gcm.data(var.name='clt',gcm=gcm,gcm.file=clt.file,lon.c,lat.c)

  rv <- list(rhs=rhs.data,
             rsds=rsds.data,
             clt=clt.data)

if (1==0) {
  gcm.dir <- paste('/storage/data/climate/downscale/CMIP5/building_code/',gcm,sep='')

  psl.file <- list.files(path=gcm.dir,pattern=paste0('psl_day_',gcm),full.name=TRUE)
  psl.data <- get.gcm.data(var.name='psl',gcm=gcm,gcm.file=psl.file,lon.c,lat.c)

  huss.file <- list.files(path=gcm.dir,pattern=paste0('huss_day_',gcm),full.name=TRUE)
  huss.data <- get.gcm.data(var.name='huss',gcm=gcm,gcm.file=huss.file,lon.c,lat.c)

  uas.file <- list.files(path=gcm.dir,pattern=paste0('uas_day_',gcm),full.name=TRUE)
  uas.data <- get.gcm.data(var.name='uas',gcm=gcm,gcm.file=uas.file,lon.c,lat.c)

  vas.file <- list.files(path=gcm.dir,pattern=paste0('vas_day_',gcm),full.name=TRUE)
  vas.data <- get.gcm.data(var.name='vas',gcm=gcm,gcm.file=vas.file,lon.c,lat.c)

  rv <- list(psl=psl.data,
             huss=huss.data,
             uas=uas.data,
             vas=vas.data)
}
if (1==0) {
  gcm.dir <- paste('/storage/data/climate/downscale/BCCAQ2/CMIP5/',gcm,sep='')

  pr.file <- list.files(path=gcm.dir,pattern=paste0('pr_day_',gcm),full.name=TRUE)
  pr.data <- get.gcm.data(var.name='pr',gcm=gcm,gcm.file=pr.file,lon.c,lat.c)

  tx.file <- list.files(path=gcm.dir,pattern=paste0('tasmax_day_',gcm),full.name=TRUE)
  tx.data <- get.gcm.data(var.name='tasmax',gcm=gcm,gcm.file=tx.file,lon.c,lat.c)

  tn.file <- list.files(path=gcm.dir,pattern=paste0('tasmin_day_',gcm),full.name=TRUE)
  tn.data <- get.gcm.data(var.name='tasmin',gcm=gcm,gcm.file=tn.file,lon.c,lat.c)

  rv <- list(pr=pr.data,
             tasmax=tx.data,
             tasmin=tn.data)
}
  return(rv)
         
}


##BCCAQ Data
gather.bccaq.data <- function(gcm,lon.c,lat.c,scenario) {

  bccaq.dir <- paste('/storage/data/scratch/ssobie/bccaq_gcm_bc_subset/',gcm,sep='')

  tasmax.data <- get.bccaq.data('tasmax',gcm,scenario,bccaq.dir,lon.c,lat.c)
  tasmin.data <- get.bccaq.data('tasmin',gcm,scenario,bccaq.dir,lon.c,lat.c)
  tas.data <- tasmax.data
  tas.data$data <- (tasmax.data$data+tasmin.data$data)/2
  pr.data <- get.bccaq.data('pr',gcm,scenario,bccaq.dir,lon.c,lat.c)

  rv <- list(tasmax=tasmax.data,        
             tasmin=tasmin.data,
             tas=tas.data,
             pr=pr.data)
  return(rv)
}


##Snow Model Data
gather.snow.data <- function(gcm,lon.c,lat.c,scenario) {

  snow.dir <- paste('/storage/data/projects/rci/data/assessments/high_res_data/bccaq_gcm_nanaimo_subset/rcp85/snow/',gcm,sep='')

  snow.data <- get.snow.data('snowdepth',gcm,scenario,snow.dir,lon.c,lat.c)

  rv <- list(snow=snow.data)
  return(rv)
}

##************************************************************************
##************************************************************************
##Load GCM BasedData


##gcm <- 'ACCESS1-0'
##scenario <- 'rcp45'

##lon.c <- -123.969357
##lat.c <- 49.184737

##BCCAQ Data
##bccaq.data <- gather.bccaq.data(gcm,lon.c,lat.c,scenario)

##gcm.list <- c('ACCESS1-0',
##              'CanESM2',
##              'CCSM4',
##              'CNRM-CM5',
##              'CSIRO-Mk3-6-0',
##              'GFDL-ESM2G',
##              'HadGEM2-CC',
##              'HadGEM2-ES',
##              'inmcm4',
##              'MIROC5',
##              'MPI-ESM-LR',
##              'MRI-CGCM3')


##GCM other variables
##scenario <- 'rcp85'

##save.dir <- '/storage/data/climate/downscale/CMIP5/building_code/data_files/'
##for (gcm in gcm.list) {
##    gcm.file <- paste0(save.dir,'nanaimo_hospital_gcm_',gcm,'_',scenario,'_pr_tx_tn.RData') 
##    gcm.data <- gather.gcm.data(gcm,lon.c,lat.c,scenario) 
##    save(gcm.data,file=gcm.file)
##}
