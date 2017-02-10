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
  data.subset <- data[st:en]
  dates.subset <- dates[st:en]
  rv <- data.subset
  return(rv)
}

get.gcm.data <- function(var.name,gcm,gcm.file,cell.subset) {

  print(gcm.file)
  nc <- nc_open(gcm.file)
  gcm.data <- ncvar_get(nc,start=c(cell.subset,1),count=c(1,1,-1))
  gcm.time <- get.time.series(nc)
  gcm.subset <- read.data.subset(gcm.data,gcm.time,'1951-2100')
  time.subset <- as.Date(read.data.subset(gcm.time,gcm.time,'1951-2100'))
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

get.bccaq.data <- function(var.name,gcm,scenario,data.dir,cell.subset) {

   gcm.files <- list.files(path=data.dir,
                           pattern=paste(var.name,'_day_BCCAQ_',gcm,'_',scenario,sep=''),full.name=TRUE)
   print(gcm.files)
   past.file <- gcm.files[grep('1951-2000',gcm.files)]
   proj.file <- gcm.files[grep('2001-2100',gcm.files)]
   nc.past <- nc_open(past.file)
   nc.proj <- nc_open(proj.file)
   past.subset <- ncvar_get(nc.past,start=c(cell.subset,1),count=c(1,1,-1))
   past.time <- get.time.series(nc.past)
   proj.subset <- ncvar_get(nc.proj,start=c(cell.subset,1),count=c(1,1,-1))
   proj.time <- get.time.series(nc.proj)
   comb.subset <- c(past.subset,proj.subset)
   time.subset <- as.Date(c(past.time,proj.time))

   nc_close(nc.past)
   nc_close(nc.proj)

   rv <- list(data=comb.subset,time=time.subset)
   return(rv)
}



##************************************************************************
##************************************************************************
##Load GCM BasedData

cell <- c(85,50)
gcm <- 'CanESM2'
scenario <- 'rcp85'

gcm.dir <- paste('/storage/data/climate/downscale/BCCAQ2/CMIP5/',gcm,sep='')

psl.file <- paste(gcm.dir,'/psl_day_',gcm,'_historical+rcp85_r1i1p1_18500101-21001231.nc',sep='')
psl.data <- get.gcm.data(var.name='psl',gcm=gcm,gcm.file=psl.file,cell)

huss.file <- paste(gcm.dir,'/huss_day_',gcm,'_historical+rcp85_r1i1p1_18500101-21001231.nc',sep='')
huss.data <- get.gcm.data(var.name='huss',gcm=gcm,gcm.file=huss.file,cell)

uas.file <- paste(gcm.dir,'/uas_day_',gcm,'_historical+rcp85_r1i1p1_18500101-21001231.nc',sep='')
uas.data <- get.gcm.data(var.name='uas',gcm=gcm,gcm.file=uas.file,cell)

vas.file <- paste(gcm.dir,'/vas_day_',gcm,'_historical+rcp85_r1i1p1_18500101-21001231.nc',sep='')
vas.data <- get.gcm.data(var.name='vas',gcm=gcm,gcm.file=vas.file,cell)

##BCCAQ Data
cell <- c(193,15)
bccaq.dir <- paste('/storage/data/scratch/ssobie/bccaq_gcm_bc_subset/',gcm,sep='')

tasmax.data <- get.bccaq.data('tasmax',gcm,scenario,bccaq.dir,cell)
tasmin.data <- get.bccaq.data('tasmin',gcm,scenario,bccaq.dir,cell)
tas.data <- tasmax.data
tas.data$data <- (tasmax.data$data+tasmin.data$data)/2
pr.data <- get.bccaq.data('pr',gcm,scenario,bccaq.dir,cell)



