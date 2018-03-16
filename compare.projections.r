##Script to compare projected changes of TAS and Precip between
##the CanRCM4 cell and the corresponding 10km cells

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

get.data.subset <- function(data,dates,interval) {

  bnds <- strsplit(interval,'-')[[1]]
  st <- head(grep(bnds[1],dates),1)
  en <- tail(grep(bnds[2],dates),1)
  data.subset <- data[st:en]
  dates.subset <- dates[st:en]
  rv <- data.subset
##  rv <- list(data=data.subset,
##             dates=dates.subset)
  return(rv)
}

fill.series <- function(data,dates,template) {
  rv <- rep(NA,length(template))
  rv[template %in% dates] <- data[template %in% dates]            
  return(rv)
}

get.gcm.data <- function(var.name,gcm.list,data.dir) {

  lon.i <- -123.968641
  lat.i <- 49.185618

  avg.subset <- vector(length=length(gcm.list),mode='list')
  time.subset <- vector(length=length(gcm.list),mode='list')
     for (i in seq_along(gcm.list)) {
         gcm <- gcm.list[i]
         gcm.files <- list.files(path=paste(data.dir,gcm,sep=''),
                                 pattern=paste(var.name,'_day_',sep=''),full.name=TRUE)
         print(gcm.files)                               
         past.file <- gcm.files[grep('rcp85',gcm.files)]

         nc.past <- nc_open(past.file)       
         lon <- ncvar_get(nc.past,'lon')
         lon <- ((lon + 180) %% 360) - 180
         lat <- ncvar_get(nc.past,'lat')
         lon.ix <- which.min(abs(lon.i-lon))
         lat.ix <- which.min(abs(lat.i-lat))

         past.subset <- ncvar_get(nc.past,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
         if (grepl('tas',var.name)) {
            past.subset <- past.subset-273.15
         }
         if (grepl('pr',var.name)) {
            past.subset <- past.subset*86400
         }

         past.time <- get.time.series(nc.past)
         avg.subset[[i]] <- past.subset
         time.subset[[i]] <- past.time
         nc_close(nc.past)
   }

  rv <- list(data=avg.subset,time=time.subset)
  return(rv)
}


get.bccaq.data <- function(var.name,gcm.list,data.dir,cell.subset) {

  avg.subset <- vector(length=length(gcm.list),mode='list')
  time.subset <- vector(length=length(gcm.list),mode='list')
     for (i in seq_along(gcm.list)) {
         gcm <- gcm.list[i]
         gcm.files <- list.files(path=paste(data.dir,gcm,sep=''),
                                 pattern=paste(var.name,'_day_BCCAQ',sep=''),full.name=TRUE)
         print(gcm.files)                               
         past.file <- gcm.files[grep('1951-2000',gcm.files)]
         proj.file <- gcm.files[grep('2001-2100',gcm.files)]
         nc.past <- nc_open(past.file)       
         nc.proj <- nc_open(proj.file)       
         past.subset <- lapply(cell.subset,function(x,nc){return(ncvar_get(nc,start=c(x,1),count=c(1,1,-1)))},nc.past)
         past.time <- get.time.series(nc.past)
         proj.subset <- lapply(cell.subset,function(x,nc){return(ncvar_get(nc,start=c(x,1),count=c(1,1,-1)))},nc.proj)
         proj.time <- get.time.series(nc.proj)
         comb.subset <- mapply(c,past.subset,proj.subset)
         avg.subset[[i]] <- apply(comb.subset,1,mean,na.rm=T)
         time.subset[[i]] <- c(past.time,proj.time)

         nc_close(nc.past)
         nc_close(nc.proj)
   }

  rv <- list(data=avg.subset,time=time.subset)
  return(rv)
}

get.season.data <- function(data,dates,fxn) {
  annual.fac <- as.factor(format(as.Date(dates),'%Y'))
  monthly.fac <- as.factor(format(as.Date(dates),'%m'))
  seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
  seasonal.fac <- factor(seasons[monthly.fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))
  ##seas.avg <- tapply(data,seasonal.fac,fxn,na.rm=TRUE)
  seas.yr.avg <- tapply(data,list(annual.fac,seasonal.fac),fxn,na.rm=TRUE)
  seas.avg <- apply(seas.yr.avg,2,mean,na.rm=T)

  return(seas.avg)

}

get.rcm.data <- function(var.name,cell,read.dir) {

   load(paste(read.dir,var.name,'_NANAIMO_',cell[1],'_',cell[2],'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
   data.rcm <- data

   if (var.name=='tas')
      data.rcm <- data.rcm - 273
   if (var.name=='pr')
      data.rcm[data.rcm < 0.01] <- 0
 
   return(data.rcm)              
}

get.rcm.seasonal <- function(var.name,data.rcm,rcm.dates) {
 
   dates.1980s <- get.data.subset(rcm.dates,rcm.dates,'1971-2000')
   dates.2050s <- get.data.subset(rcm.dates,rcm.dates,'2041-2070')

   rcm.1980s <- get.data.subset(data.rcm,rcm.dates,'1971-2000')
   rcm.2020s <- get.data.subset(data.rcm,rcm.dates,'2011-2040')
   rcm.2050s <- get.data.subset(data.rcm,rcm.dates,'2041-2070')
   rcm.2080s <- get.data.subset(data.rcm,rcm.dates,'2071-2100')

   if (var.name=='pr') {
     rcm.anoms <- (mean(rcm.2050s)-mean(rcm.1980s))/mean(rcm.1980s)*100  
     rcm.seas.1980s <- get.season.data(rcm.1980s,dates.1980s,'sum')
     rcm.seas.2050s <- get.season.data(rcm.2050s,dates.2050s,'sum')
     rcm.seas.anoms <- (rcm.seas.2050s - rcm.seas.1980s)/rcm.seas.1980s*100
     rcm.seasmax.1980s <- get.season.data(rcm.1980s,dates.1980s,'max')
     rcm.seasmax.2050s <- get.season.data(rcm.2050s,dates.2050s,'max')
     rcm.seasmax.anoms <- (rcm.seasmax.2050s - rcm.seasmax.1980s)/rcm.seasmax.1980s*100
   } else {
     rcm.anoms <- mean(rcm.2050s)-mean(rcm.1980s)
     rcm.seas.1980s <- get.season.data(rcm.1980s,dates.1980s,'mean')
     rcm.seas.2050s <- get.season.data(rcm.2050s,dates.2050s,'mean')   
     rcm.seas.anoms <- rcm.seas.2050s - rcm.seas.1980s        
     rcm.seasmax.1980s <- get.season.data(rcm.1980s,dates.1980s,'max')
     rcm.seasmax.2050s <- get.season.data(rcm.2050s,dates.2050s,'max')
     rcm.seasmax.anoms <- rcm.seasmax.2050s - rcm.seasmax.1980s
   }

   rv <- round(rbind(c(rcm.seas.1980s,rcm.seas.2050s,rcm.seas.anoms),
               c(rcm.seasmax.1980s,rcm.seasmax.2050s,rcm.seasmax.anoms)),1)
   return(rv)              
}

gcm.or.rcm.data <- function(pr.data,tasmax.data,tasmin.data,tas.data) {


   bccaq.1980s <- mapply(get.data.subset,tasmax.data$time,tasmax.data$time,MoreArgs=list('1971-2000'))
   bccaq.2050s <- mapply(get.data.subset,tasmax.data$time,tasmax.data$time,MoreArgs=list('2041-2070'))

   tas.1980s <- mapply(get.data.subset,tas.data,tasmax.data$time,MoreArgs=list('1971-2000'))
   tas.2020s <- mapply(get.data.subset,tas.data,tasmax.data$time,MoreArgs=list('2011-2040'))
   tas.2050s <- mapply(get.data.subset,tas.data,tasmax.data$time,MoreArgs=list('2041-2070'))
   tas.2080s <- mapply(get.data.subset,tas.data,tasmax.data$time,MoreArgs=list('2071-2099'))

   tasmax.1980s <- mapply(get.data.subset,tasmax.data$data,tasmax.data$time,MoreArgs=list('1971-2000'))
   tasmax.2020s <- mapply(get.data.subset,tasmax.data$data,tasmax.data$time,MoreArgs=list('2011-2040'))
   tasmax.2050s <- mapply(get.data.subset,tasmax.data$data,tasmax.data$time,MoreArgs=list('2041-2070'))
   tasmax.2080s <- mapply(get.data.subset,tasmax.data$data,tasmax.data$time,MoreArgs=list('2071-2099'))

   pr.1980s <- mapply(get.data.subset,pr.data$data,tasmax.data$time,MoreArgs=list('1971-2000'))
   pr.2020s <- mapply(get.data.subset,pr.data$data,tasmax.data$time,MoreArgs=list('2011-2040'))
   pr.2050s <- mapply(get.data.subset,pr.data$data,tasmax.data$time,MoreArgs=list('2041-2070'))
   pr.2080s <- mapply(get.data.subset,pr.data$data,tasmax.data$time,MoreArgs=list('2071-2099'))

   tas.clim.1980s <- unlist(lapply(tas.1980s,mean,na.rm=T))
   tas.clim.2050s <- unlist(lapply(tas.2050s,mean,na.rm=T))
   tas.anoms <- tas.clim.2050s - tas.clim.1980s

   pr.clim.1980s <- unlist(lapply(pr.1980s,function(x){mean(x[x>0],na.rm=T)}))
   pr.clim.2050s <- unlist(lapply(pr.2050s,function(x){mean(x[x>0],na.rm=T)}))
   pr.anoms <- (pr.clim.2050s - pr.clim.1980s)/pr.clim.1980s*100

   tas.seas.1980s <- mapply(get.season.data,tas.1980s,bccaq.1980s,'mean')
   tas.seas.2050s <- mapply(get.season.data,tas.2050s,bccaq.2050s,'mean')
   pr.seas.1980s <- mapply(get.season.data,pr.1980s,bccaq.1980s,'sum')
   pr.seas.2050s <- mapply(get.season.data,pr.2050s,bccaq.2050s,'sum')

   tasmax.seasmax.1980s <- mapply(get.season.data,tasmax.1980s,bccaq.1980s,'max')
   tasmax.seasmax.2050s <- mapply(get.season.data,tasmax.2050s,bccaq.2050s,'max')

   pr.seasmax.1980s <- mapply(get.season.data,pr.1980s,bccaq.1980s,'max')
   pr.seasmax.2050s <- mapply(get.season.data,pr.2050s,bccaq.2050s,'max')

   tas.seas.anoms <- tas.seas.2050s - tas.seas.1980s
   pr.seas.anoms <- (pr.seas.2050s - pr.seas.1980s)/pr.seas.1980s*100

   tasmax.seasmax.anoms <- tasmax.seasmax.2050s - tasmax.seasmax.1980s
   pr.seasmax.anoms <- (pr.seasmax.2050s - pr.seasmax.1980s)/pr.seasmax.1980s*100


   rv <- list(tas.seas.1980s=tas.seas.1980s,tas.seas.2050s=tas.seas.2050s,
              pr.seas.1980s=pr.seas.1980s,pr.seas.2050s=pr.seas.2050s,
              tasmax.seasmax.1980s=tasmax.seasmax.1980s,tasmax.seasmax.2050s=tasmax.seasmax.2050s,
              pr.seasmax.1980s=pr.seasmax.1980s,pr.seasmax.2050s=pr.seasmax.2050s,
              tas.seas.anoms=tas.seas.anoms,pr.seas.anoms=pr.seas.anoms,
              tasmax.seasmax.anoms=tasmax.seasmax.anoms,pr.seasmax.anoms=pr.seasmax.anoms)
   return(rv)
}


##************************************************************************
##************************************************************************
##Load RCM Data
read.dir <- '/storage/data/climate/downscale/BCCAQ2/CanRCM4/'

rcm.dates <- seq(from=as.Date('1950-01-01'),by='day',to=as.Date('2100-12-31'))
flag <- grep('*-02-29',rcm.dates)
rcm.dates <- rcm.dates[-flag]

cell.list <- list(c(77,153),
                  c(76,153),    
                  c(76,152),   
                  c(75,152),
                  c(75,151),
                  c(76,151),
                  c(76,150))

cell.list <- list(c(77,153),
                  c(76,150))


rlen <- length(cell.list)
tas.list <- vector(length=rlen,mode='list')
tasmax.list <- vector(length=rlen,mode='list')
tasmin.list <- vector(length=rlen,mode='list')
pr.list <- vector(length=rlen,mode='list')

tas.rcm.series <- vector(length=rlen,mode='list')
tasmax.rcm.series <- vector(length=rlen,mode='list')
tasmin.rcm.series <- vector(length=rlen,mode='list')
pr.rcm.series <- vector(length=rlen,mode='list')

for (i in 1:rlen) {
    cell <- cell.list[[i]]    
    pr.rcm.series[[i]] <- get.rcm.data('pr',cell,read.dir)
    pr.list[[i]] <- get.rcm.seasonal('pr',pr.rcm.series[[i]],rcm.dates)
    tas.rcm.series[[i]] <- get.rcm.data('tas',cell,read.dir)
    tas.list[[i]] <- get.rcm.seasonal('tas',tas.rcm.series[[i]],rcm.dates)
    tasmax.rcm.series[[i]] <- get.rcm.data('tasmax',cell,read.dir)
    tasmax.list[[i]] <- get.rcm.seasonal('tasmax',tasmax.rcm.series[[i]],rcm.dates)
    tasmin.rcm.series[[i]] <- get.rcm.data('tasmin',cell,read.dir)
    tasmin.list[[i]] <- get.rcm.seasonal('tasmin',tasmin.rcm.series[[i]],rcm.dates)

}

##---------------------------------------------------------------------
##Load BCCAQ Data
data.dir <- '/storage/data/scratch/ssobie/bccaq_gcm_nanaimo_subset/'
gcm.dir <- '/storage/data/climate/downscale/BCCAQ2/CMIP5/'
gcm.list <- c('ACCESS1-0',
              'CanESM2',
              'CCSM4',
              'CNRM-CM5',
              'CSIRO-Mk3-6-0',
              'GFDL-ESM2G',
              'HadGEM2-CC',
              'HadGEM2-ES',
              'inmcm4',
              'MIROC5',
              'MPI-ESM-LR',
              'MRI-CGCM3')

cell.subset <- list(c(3,2),c(3,3),c(3,4),c(3,5),
                    c(4,2),c(4,3),c(4,4),
                    c(5,2),c(5,3))

if (1==0) {
pr.data <- get.bccaq.data('pr',gcm.list,data.dir,cell.subset)
tasmax.data <- get.bccaq.data('tasmax',gcm.list,data.dir,cell.subset)
tasmin.data <- get.bccaq.data('tasmin',gcm.list,data.dir,cell.subset)
tas.data <- mapply(FUN=function(x,y){(x+y)/2},tasmax.data$data,tasmin.data$data)

bccaq.vars <- gcm.or.rcm.data(pr.data,tasmax.data,tasmin.data,tas.data)

##GCM Data
gcm.pr.data <- get.gcm.data('pr',gcm.list,gcm.dir)
gcm.tasmax.data <- get.gcm.data('tasmax',gcm.list,gcm.dir)
gcm.tasmin.data <- get.gcm.data('tasmin',gcm.list,gcm.dir)
gcm.tas.data <- mapply(FUN=function(x,y){(x+y)/2},gcm.tasmax.data$data,gcm.tasmin.data$data)

gcm.vars <- gcm.or.rcm.data(gcm.pr.data,gcm.tasmax.data,gcm.tasmin.data,gcm.tas.data)

}

##**********************************************************************************
##**********************************************************************************
##Plotting

plot.dir <- '/storage/data/projects/rci/building_code/plots/'

if (1==0) {
##Climatologies
plot.file <- paste(plot.dir,'gcm.rcm.bccaq.mean.tas.pr.1980s.climatologies_scatter.png',sep='')
seas.title <- c('Winter Mean','Spring Mean','Summer Mean','Fall Mean')
png(file=plot.file,width=900,height=900)
j <- 2
xlims <- list(c(-5,5.5),c(2,9.5),c(12,19),c(4,11))
ylims <- list(c(500,900),c(250,600),c(40,260),c(300,900))
par(mfrow=c(2,2),mar=c(5,5,2,2))
for (i in 1:4) {
  plot(bccaq.vars$tas.seas.1980s[i,],bccaq.vars$pr.seas.1980s[i,],
        pch=15,col='red',xlab='Annual Avg TAS (degC)',ylab='Annual Total PR (mm)',main=seas.title[i],
        xlim=xlims[[i]], ##c(min(tas.seas.1980s[i,])-1,max(tas.seas.1980s[i,])+1),
        ylim=ylims[[i]], ##c(min(pr.seas.1980s[i,])-200,max(pr.seas.1980s[i,])+200),
        cex=3,cex.main=2,cex.lab=2,cex.axis=2)
  points(tas.list[[j]][1,i],pr.list[[j]][1,i],pch=17,col='blue',cex=3)
  points(gcm.vars$tas.seas.1980s[i,],gcm.vars$pr.seas.1980s[i,],col='green',cex=3,pch=16)
  if (i==4) {
  legend('topright',legend=c('GCM ~150km','BCCAQ 30km','RCM 25km'),col=c('green','red','blue'),pch=c(16,15,17),cex=2)
  }
}
dev.off()




plot.file <- paste(plot.dir,'gcm.rcm.max.tas.pr.1980s.climatologies_scatter.png',sep='')
seas.title <- c('Winter Max','Spring Max','Summer Max','Fall Max')
xlims <- list(c(1.5,13.5),c(17,28),c(26,38),c(22,31))
ylims <- list(c(30,70),c(20,50),c(10,30),c(25,70))
png(file=plot.file,width=900,height=900)
j<-2
par(mfrow=c(2,2),mar=c(5,5,2,2))
for (i in 1:4) {
  plot(bccaq.vars$tasmax.seasmax.1980s[i,],bccaq.vars$pr.seasmax.1980s[i,],
        pch=15,col='red',xlab='Annual Avg TASMAX (degC)',ylab='Annual Avg Max Precip (mm)',main=seas.title[i],
        xlim=xlims[[i]], ##c(min(tasmax.seasmax.1980s[i,])-1,max(tasmax.seasmax.1980s[i,])+1),
        ylim=ylims[[i]], ##c(min(pr.seasmax.1980s[i,])-40,max(pr.seasmax.1980s[i,])+40),
        cex=3,cex.main=2,cex.lab=2,cex.axis=2)
  points(tasmax.list[[j]][2,i],pr.list[[j]][2,i],pch=17,col='blue',cex=3)
  points(gcm.vars$tasmax.seasmax.1980s[i,],gcm.vars$pr.seasmax.1980s[i,],col='green',cex=3,pch=16)
  if (i==4) {
  legend('bottomright',legend=c('GCM ~150km','BCCAQ 30km','RCM 25km'),col=c('green','red','blue'),pch=c(16,15,17),cex=2)
  }
}
dev.off()
}


if (1==0) {

plot.file <- paste(plot.dir,'gcm.rcm.max.tas.pr.1980s.anomalies_scatter.png',sep='')
seas.title <- c('Winter Max Change','Spring Max Change','Summer Max Change','Fall Max Change')
xlims <- list(c(-1,7.5),c(-0.5,6),c(0,7),c(0,5))
ylims <- list(c(-10,40),c(-2,60),c(-40,30),c(-20,60))
j<-2
png(file=plot.file,width=900,height=900)
par(mfrow=c(2,2),mar=c(5,5,2,2))
for (i in 1:4) {
  plot(bccaq.vars$tasmax.seasmax.anoms[i,],bccaq.vars$pr.seasmax.anoms[i,],
        pch=15,col='red',xlab='Annual Avg TASMAX (degC)',ylab='Annual Average MAX PR (%)',main=seas.title[i],
        xlim=xlims[[i]], ##c(min(tasmax.seasmax.anoms[i,])-1,max(tasmax.seasmax.anoms[i,])+1),
        ylim=ylims[[i]], ##c(min(pr.seasmax.anoms[i,])-30,max(pr.seasmax.anoms[i,])+30),
        cex=3,cex.main=2,cex.lab=2,cex.axis=2)
    points(tasmax.list[[j]][1,i+8],pr.list[[j]][1,i+8],pch=17,col='blue',cex=3)
    points(gcm.vars$tasmax.seasmax.anoms[i,],gcm.vars$pr.seasmax.anoms[i,],col='green',cex=3,pch=16)
    if (i==4) {
      legend('topleft',legend=c('GCM ~150km','BCCAQ 30km','RCM 25km'),col=c('green','red','blue'),pch=c(16,15,17),cex=2)
    }
}
dev.off()
}

if (1==1) {
##Projected Changes
plot.file <- paste(plot.dir,'gcm.rcm.mean.tas.pr.2050s.anoms_scatter.png',sep='')
seas.title <- c('Winter Change','Spring Change','Summer Change','Fall Change')
xlims <- list(c(1,4),c(1,5.5),c(1,6),c(0,5))
ylims <- list(c(-10,20),c(-20,20),c(-60,40),c(-20,40))
j<-2
png(file=plot.file,width=900,height=900)
par(mfrow=c(2,2),mar=c(5,5,2,2))
for (i in 1:4) {
  plot(bccaq.vars$tas.seas.anoms[i,],bccaq.vars$pr.seas.anoms[i,],
        pch=15,col='red',xlab='Annual Avg TAS (degC)',ylab='Annual Total PR (%)',main=seas.title[i],
        xlim=xlims[[i]], ##c(min(tas.seas.anoms[i,])-1,max(tas.seas.anoms[i,])+1),
        ylim=ylims[[i]], ##c(min(pr.seas.anoms[i,])-30,max(pr.seas.anoms[i,])+30),
        cex=3,cex.main=2,cex.lab=2,cex.axis=2)
    points(tas.list[[j]][2,i+8],pr.list[[j]][2,i+8],pch=17,col='blue',cex=3)
    points(gcm.vars$tas.seas.anoms[i,],gcm.vars$pr.seas.anoms[i,],col='green',cex=3,pch=16)
    if (i==4) {
      legend('topleft',legend=c('GCM ~150km','BCCAQ 30km','RCM 25km'),col=c('green','red','blue'),pch=c(16,15,17),cex=2)
    }
}
dev.off()
browser()
}



if (1==0) {
##Climatologies
plot.file <- paste(plot.dir,'mean.tas.pr.1980s.climatologies_',cell.name,'_comb.png',sep='')
seas.title <- c('Winter Mean','Spring Mean','Summer Mean','Fall Mean')
png(file=plot.file,width=900,height=900)
par(mfrow=c(2,2),mar=c(5,5,2,2))
for (i in 1:4) {
  plot(tas.seas.1980s[i,],pr.seas.1980s[i,],pch=15,col='red',xlab='Annual Avg TAS (degC)',ylab='Annual Total PR (mm)',main=seas.title[i],
        xlim=c(min(tas.seas.1980s[i,])-4,max(tas.seas.1980s[i,])+4),
        ylim=c(min(pr.seas.1980s[i,])-600,max(pr.seas.1980s[i,])+600),
        cex=2,cex.main=2,cex.lab=2,cex.axis=2)
  for (j in 1:rlen) {      
    points(tas.list[[j]][1,i],pr.list[[j]][1,i],pch=15,col='blue',cex=2)
    text(x=tas.list[[j]][1,i],y=pr.list[[j]][1,i],labels=paste(cell.list[[j]],collapse=','))    
  }
}
dev.off()


plot.file <- paste(plot.dir,'max.tas.pr.1980s.climatologies_',cell.name,'_comb.png',sep='')
seas.title <- c('Winter Max','Spring Max','Summer Max','Fall Max')
png(file=plot.file,width=900,height=900)
par(mfrow=c(2,2),mar=c(5,5,2,2))
for (i in 1:4) {
  plot(tasmax.seasmax.1980s[i,],pr.seasmax.1980s[i,],pch=15,col='red',xlab='Annual Avg TAS (degC)',ylab='Annual Total PR (mm)',main=seas.title[i],
        xlim=c(min(tasmax.seasmax.1980s[i,])-1,max(tasmax.seasmax.1980s[i,])+1),
        ylim=c(min(pr.seasmax.1980s[i,])-40,max(pr.seasmax.1980s[i,])+40),
        cex=2,cex.main=2,cex.lab=2,cex.axis=2)
    points(tasmax.list[[j]][2,i],pr.list[[j]][2,i],pch=15,col='blue',cex=2)
}
dev.off()

plot.file <- paste(plot.dir,'max.tas.pr.1980s.anomalies_',cell.name,'.png',sep='')
seas.title <- c('Winter Max Change','Spring Max Change','Summer Max Change','Fall Max Change')
png(file=plot.file,width=900,height=900)
par(mfrow=c(2,2),mar=c(5,5,2,2))
for (i in 1:4) {
  plot(tasmax.seasmax.anoms[i,],pr.seasmax.anoms[i,],pch=15,col='red',xlab='Annual Avg TASMAX (degC)',ylab='Annual Average MAX PR (%)',main=seas.title[i],
        xlim=c(min(tasmax.seasmax.anoms[i,])-1,max(tasmax.seasmax.anoms[i,])+1),
        ylim=c(min(pr.seasmax.anoms[i,])-30,max(pr.seasmax.anoms[i,])+30),
        cex=2,cex.main=2,cex.lab=2,cex.axis=2)
  for (j in 1:rlen) {      
    points(tasmax.list[[j]][1,i+8],pr.list[[j]][1,i+8],pch=15,col='blue',cex=2)
    text(x=tasmax.list[[j]][1,i+8],y=pr.list[[j]][1,i+8],labels=paste(cell.list[[j]],collapse=','))    
  }
}
dev.off()


##Projected Changes
plot.file <- paste(plot.dir,'mean.tas.pr.2050s.anoms_',cell.name,'.png',sep='')
seas.title <- c('Winter Change','Spring Change','Summer Change','Fall Change')
png(file=plot.file,width=900,height=900)
par(mfrow=c(2,2),mar=c(5,5,2,2))
for (i in 1:4) {
  plot(tas.seas.anoms[i,],pr.seas.anoms[i,],pch=15,col='red',xlab='Annual Avg TAS (degC)',ylab='Annual Total PR (%)',main=seas.title[i],
        xlim=c(min(tas.seas.anoms[i,])-1,max(tas.seas.anoms[i,])+1),
        ylim=c(min(pr.seas.anoms[i,])-30,max(pr.seas.anoms[i,])+30),
        cex=2,cex.main=2,cex.lab=2,cex.axis=2)
  for (j in 1:rlen) {      
    points(tas.list[[j]][2,i+8],pr.list[[j]][2,i+8],pch=15,col='blue',cex=2)
    text(x=tas.list[[j]][2,i+8],y=pr.list[[j]][2,i+8],labels=paste(cell.list[[j]],collapse=','))    
  }
}
dev.off()

}

if (1==0) {
date.max <- bccaq.1980s[[1]]
matrix.1980s <- mapply(fill.series,tas.1980s,bccaq.1980s,MoreArgs=list(date.max))
mean.1980s <- apply(matrix.1980s,1,mean,na.rm=T)

plot(tas.rcm.series[[7]][['date1980s']],tas.rcm.series[[7]][['hist']],ylim=c(-20,35))
##lapply(tas.1980s,lines)
##lines(tas.rcm.series[[7]][['date1980s']],tas.rcm.series[[7]][['hist']],lwd=3,col='green')
lines(as.Date(date.max),mean.1980s,lwd=3,col='red')

mon.tas <- function(x,y) {
   rv <- tapply(x,as.factor(format(as.Date(y),'%Y')),mean)
   return(rv)
} 
mon.pr <- function(x,y) {
   rv <- tapply(x,as.factor(format(as.Date(y),'%Y')),sum)
   return(rv)
} 

tas.mon <- mapply(mon.tas,tas.1980s,bccaq.1980s)
pr.mon <- mapply(mon.pr,pr.data$data,pr.data$time)

pr.rcm.mon <- mapply(mon.pr,pr.rcm.1980s,rcm.1980s)

}

##***********************************************************
##***********************************************************
if (1==0) {

ix <- 2
ann.mean <- c()
yr.rcm.fac <- as.factor(format(rcm.dates,'%Y'))
rcm.ann.sum <- tapply(pr.rcm.series[[ix]],yr.rcm.fac,sum)

plot.file <- paste(plot.dir,'pr.series.total.png',sep='')
png(file=plot.file,width=900,height=600)
par(mar=c(5,5,2,2))
plot(c(),col='white',xlab='Years',ylab='Precipitation (mm)',main='Annual Total Precipitation',
     xlim=c(as.Date('1950-01-01'),as.Date('2100-01-01')),
     ylim=c(1000,3500),
     cex=2,cex.main=2,cex.lab=2,cex.axis=2,axes=F) 
axis(1,at=as.Date(paste(seq(1950,2100,by=25),'-01-01',sep='')),labels=seq(1950,2100,by=25),cex.axis=2,cex.lab=2)
axis(2,at=seq(1000,3500,by=500),label=seq(1000,3500,by=500),cex.axis=2,cex.lab=2)
for (i in 1:12) {
  yr.fac <- as.factor(format(as.Date(pr.data$time[[i]]),'%Y'))
  ann.sum <- tapply(pr.data$data[[i]],yr.fac,sum,na.rm=T)
  print(length(ann.sum))
  ann.mean <- rbind(ann.mean,ann.sum)
  print(dim(ann.mean))
  lines(as.Date(paste(levels(yr.fac),'-01-01',sep='')),ann.sum,col='lightgray',lwd=2)
}
ann.mean[7:8,150] <- NA
abline(h=seq(1000,3500,by=500),lty=2,col='lightgray')
lines(as.Date(paste(levels(yr.fac),'-01-01',sep='')),apply(ann.mean,2,mean),col='blue',lwd=4)
lines(as.Date(paste(levels(yr.rcm.fac),'-01-01',sep='')),rcm.ann.sum,col='red',lwd=4)
box(which='plot')
dev.off()

##***********************************************************
ann.mean <- c()
yr.rcm.fac <- as.factor(format(rcm.dates,'%Y'))
rcm.ann.sum <- tapply(pr.rcm.series[[ix]],yr.rcm.fac,max)

plot.file <- paste(plot.dir,'pr.series.max.png',sep='')
png(file=plot.file,width=900,height=600)
par(mar=c(5,5,2,2))
plot(c(),col='white',xlab='Years',ylab='Precipitation (mm)',main='Annual Max Precipitation',
     xlim=c(as.Date('1950-01-01'),as.Date('2100-01-01')),
     ylim=c(25,200),
     cex=2,cex.main=2,cex.lab=2,cex.axis=2,axes=F)
axis(1,at=as.Date(paste(seq(1950,2100,by=25),'-01-01',sep='')),labels=seq(1950,2100,by=25),cex.axis=2,cex.lab=2)
axis(2,at=seq(25,200,by=25),label=seq(25,200,by=25),cex.axis=2,cex.lab=2)
for (i in 1:12) {
  yr.fac <- as.factor(format(as.Date(pr.data$time[[i]]),'%Y'))
  ann.sum <- tapply(pr.data$data[[i]],yr.fac,max,na.rm=T)
  print(length(ann.sum))
  ann.mean <- rbind(ann.mean,ann.sum)
  print(dim(ann.mean))
  lines(as.Date(paste(levels(yr.fac),'-01-01',sep='')),ann.sum,col='lightgray',lwd=2)
}
ann.mean[7:8,150] <- NA
abline(h=seq(25,200,by=25),lty=2,col='lightgray')
lines(as.Date(paste(levels(yr.fac),'-01-01',sep='')),apply(ann.mean,2,mean),col='blue',lwd=4)
lines(as.Date(paste(levels(yr.rcm.fac),'-01-01',sep='')),rcm.ann.sum,col='red',lwd=4)
box(which='plot')
dev.off()

##***********************************************************
ann.mean <- c()
yr.rcm.fac <- as.factor(format(rcm.dates,'%Y'))
rcm.ann.mean <- tapply(tas.rcm.series[[ix]],yr.rcm.fac,mean)

plot.file <- paste(plot.dir,'tas.series.mean.png',sep='')
png(file=plot.file,width=900,height=600)
par(mar=c(5,5,2,2))
plot(c(),col='white',xlab='Years',ylab='Temperature (degC)',main='Annual Average Temperature',
     xlim=c(as.Date('1950-01-01'),as.Date('2100-01-01')),
     ylim=c(5,16),
     cex=2,cex.main=2,cex.lab=2,cex.axis=2,axes=F)
axis(1,at=as.Date(paste(seq(1950,2100,by=25),'-01-01',sep='')),labels=seq(1950,2100,by=25),cex.axis=2,cex.lab=2)
axis(2,at=seq(5,20,by=5),label=seq(5,20,by=5),cex.axis=2,cex.lab=2)
for (i in 1:12) {
  yr.fac <- as.factor(format(as.Date(tasmax.data$time[[i]]),'%Y'))
  ann.avg <- tapply(tas.data[[i]],yr.fac,mean,na.rm=T)
  ann.mean <- rbind(ann.mean,ann.avg)
  print(dim(ann.mean))
  lines(as.Date(paste(levels(yr.fac),'-01-01',sep='')),ann.avg,col='lightgray',lwd=2)
}
ann.mean[7:8,150] <- NA
abline(h=seq(5,20,by=5),lty=2,col='lightgray')
lines(as.Date(paste(levels(yr.fac),'-01-01',sep='')),apply(ann.mean,2,mean),col='orange',lwd=4)
lines(as.Date(paste(levels(yr.rcm.fac),'-01-01',sep='')),rcm.ann.mean,col='red',lwd=4)
box(which='plot')
dev.off()

##***********************************************************
ann.mean <- c()
yr.rcm.fac <- as.factor(format(rcm.dates,'%Y'))
rcm.ann.mean <- tapply(tasmax.rcm.series[[ix]],yr.rcm.fac,max)

plot.file <- paste(plot.dir,'tasmax.series.mean.png',sep='')
png(file=plot.file,width=900,height=600)
par(mar=c(5,5,2,2))
plot(c(),col='white',xlab='Years',ylab='Max Temperature (degC)',main='Annual Max Temperature',
     xlim=c(as.Date('1950-01-01'),as.Date('2100-01-01')),
     ylim=c(22,50),
     cex=2,cex.main=2,cex.lab=2,cex.axis=2,axes=F)
axis(1,at=as.Date(paste(seq(1950,2100,by=25),'-01-01',sep='')),labels=seq(1950,2100,by=25),cex.axis=2,cex.lab=2)
axis(2,at=seq(20,50,by=5),label=seq(20,50,by=5),cex.axis=2,cex.lab=2)
for (i in 1:12) {
  yr.fac <- as.factor(format(as.Date(tasmax.data$time[[i]]),'%Y'))
  ann.avg <- tapply(tasmax.data$data[[i]],yr.fac,max,na.rm=T)
  ann.mean <- rbind(ann.mean,ann.avg)
  print(dim(ann.mean))
  lines(as.Date(paste(levels(yr.fac),'-01-01',sep='')),ann.avg,col='lightgray',lwd=2)
}
ann.mean[7:8,150] <- NA
abline(h=seq(15,50,by=5),lty=2,col='lightgray')
lines(as.Date(paste(levels(yr.fac),'-01-01',sep='')),apply(ann.mean,2,mean),col='orange',lwd=4)
lines(as.Date(paste(levels(yr.rcm.fac),'-01-01',sep='')),rcm.ann.mean,col='red',lwd=4)
box(which='plot')
dev.off()



##***********************************************************
ann.mean <- c()
yr.rcm.fac <- as.factor(format(rcm.dates,'%Y'))
rcm.ann.mean <- tapply(tasmin.rcm.series[[ix]],yr.rcm.fac,min)

plot.file <- paste(plot.dir,'tasmin.series.mean.png',sep='')
png(file=plot.file,width=900,height=600)
par(mar=c(5,5,2,2))
plot(c(),col='white',xlab='Years',ylab='Min Temperature (degC)',main='Annual Min Temperature',
     xlim=c(as.Date('1950-01-01'),as.Date('2100-01-01')),
     ylim=c(-30,5),
     cex=2,cex.main=2,cex.lab=2,cex.axis=2,axes=F)
axis(1,at=as.Date(paste(seq(1950,2100,by=25),'-01-01',sep='')),labels=seq(1950,2100,by=25),cex.axis=2,cex.lab=2)
axis(2,at=seq(-35,10,by=5),label=seq(-35,10,by=5),cex.axis=2,cex.lab=2)
for (i in 1:12) {
  yr.fac <- as.factor(format(as.Date(tasmin.data$time[[i]]),'%Y'))
  ann.avg <- tapply(tasmin.data$data[[i]],yr.fac,min,na.rm=T)
  ann.mean <- rbind(ann.mean,ann.avg)
  print(dim(ann.mean))
  lines(as.Date(paste(levels(yr.fac),'-01-01',sep='')),ann.avg,col='lightgray',lwd=2)
}
ann.mean[7:8,150] <- NA
abline(h=seq(-35,10,by=5),lty=2,col='lightgray')
lines(as.Date(paste(levels(yr.fac),'-01-01',sep='')),apply(ann.mean,2,mean),col='orange',lwd=4)
lines(as.Date(paste(levels(yr.rcm.fac),'-01-01',sep='')),rcm.ann.mean,col='red',lwd=4)
box(which='plot')
dev.off()

}




