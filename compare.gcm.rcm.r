##Script to compare projected changes of TAS and Precip between
##the CanRCM4 cell and the corresponding 10km cells

library(ncdf4)
library(PCICt)
library(zyp)

get.trend <- function(avg,yrs) {

  tr.series <- zyp.trend.dataframe(as.data.frame(t(avg)), 0, "zhang", conf.intervals=T)
  ts <- yrs
  len <- length(yrs)

  series.b <- mean(avg,na.rm=T) - tr.series$trend*len/2
  series.m <- tr.series$trend
  series.sig <- tr.series$sig
  ##series.trend <- series.m*(1:len) + series.b
  series.lb <- tr.series$lbound
  series.ub <- tr.series$ubound
  rv <- list(slope=series.m,intercept=series.b,lbound=series.lb,ubound=series.ub,sig=series.sig)
}


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
  return(rv)
}

fill.series <- function(data,dates,template) {
  rv <- rep(NA,length(template))
  rv[template %in% dates] <- data[template %in% dates]            
  return(rv)
}

get.gcm.data <- function(var.name,gcm,data.dir,cell.subset) {

  gcm.file <- paste(data.dir,'/',var.name,'_day_CanESM2_historical+rcp85_r1i1p1_18500101-21001231.nc',sep='')
  print(gcm.file)
  nc <- nc_open(gcm.file)       
  gcm.data <- ncvar_get(nc,start=c(cell.subset,1),count=c(1,1,-1))
  gcm.time <- get.time.series(nc)
  gcm.subset <- get.data.subset(gcm.data,gcm.time,'1950-2100')
  time.subset <- get.data.subset(gcm.time,gcm.time,'1950-2100')
  nc_close(nc)
  if (grepl(var.name,'(tas|tasmax|tasmin)'))
     gcm.subset <- gcm.subset - 273
  if (grepl(var.name,'(pr)'))
     gcm.subset <- gcm.subset * 86400
  if (grepl(var.name,'(psl)'))
     gcm.subset <- gcm.subset / 1000
     

  rv <- list(data=gcm.subset,time=time.subset)
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
   if (var.name=='psl')
      data.rcm <- data.rcm/1000
 
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

compare.annual.series <- function(gcm.data,rcm.series,rcm.dates,plot.dir,
                                  var.name,var.units,fxn,var.range,var.seq,
                                  plot.file) {

   yr.rcm.fac <- as.factor(format(rcm.dates,'%Y'))
   yrs <- as.numeric(levels(yr.rcm.fac))
   rcm.ann.mean <- tapply(rcm.series,yr.rcm.fac,fxn)
   rcm.trend <- get.trend(rcm.ann.mean,yrs)
   rcm.trendline <- rcm.trend$slope*c(1,length(yrs)) + rcm.trend$intercept
   rcm.x <- c(yrs[1],tail(yrs,1))

   yr.fac <- as.factor(format(as.Date(gcm.data$time),'%Y'))
   ann.avg <- tapply(gcm.data$data,yr.fac,fxn,na.rm=T)
   yrs <- as.numeric(levels(yr.fac))
   gcm.trend <- get.trend(ann.avg,yrs)
   gcm.trendline <- gcm.trend$slope*c(1,length(yrs)) + gcm.trend$intercept
   gcm.x <- c(yrs[1],tail(yrs,1))

   plot.file <- paste(plot.dir,var.name,'.series.',fxn,'.png',sep='')
   png(file=plot.file,width=900,height=400)
   par(mar=c(5,5,2,2))
   plot(c(),col='white',xlab='Years',ylab=paste(toupper(var.name),' (',var.units,')',sep=''),
        main=paste('Annual ',toupper(fxn),' ',toupper(var.name),sep=''),
        xlim=c(as.Date('1950-01-01'),as.Date('2100-01-01')),
        ylim=var.range,
        cex=2,cex.main=2,cex.lab=2,cex.axis=2,axes=F)
        axis(1,at=as.Date(paste(seq(1950,2100,by=25),'-01-01',sep='')),labels=seq(1950,2100,by=25),cex.axis=2,cex.lab=2)
        axis(2,at=var.seq,label=var.seq,cex.axis=2,cex.lab=2)

        lines(as.Date(paste(levels(yr.fac),'-01-01',sep='')),ann.avg,col='blue',lwd=3)
        lines(as.Date(paste(gcm.x,'-01-01',sep='')),gcm.trendline,col='blue',lwd=3)
        abline(h=var.seq,lty=2,col='lightgray')
        lines(as.Date(paste(levels(yr.rcm.fac),'-01-01',sep='')),rcm.ann.mean,col='red',lwd=4)
        lines(as.Date(paste(rcm.x,'-01-01',sep='')),rcm.trendline,col='red',lwd=3)
        legend('topright',legend=c('GCM','RCM'),col=c('Blue','Red'),pch=15,cex=2)
        box(which='plot')
        dev.off()

   var.diff <- (ann.avg-rcm.ann.mean)/rcm.ann.mean*100
   plot.file <- paste(plot.dir,var.name,'.series.',fxn,'.difference.png',sep='')
   png(file=plot.file,width=900,height=400)
   par(mar=c(5,5,2,2))
   plot(as.Date(paste(levels(yr.fac),'-01-01',sep='')),var.diff,type='l',lwd=3,
        col='orange',xlab='Years',ylab=paste(toupper(var.name),' % Difference from RCM',sep=''),
        main=paste('Annual ',fxn,' ',toupper(var.name),' Differences',sep=''),
        cex=2,cex.main=2,cex.lab=2,cex.axis=2)
   abline(h=0,col='lightgray',lty=2)
   box(which='plot')
   dev.off()

}


##************************************************************************
##************************************************************************
prct <- FALSE

var.name <- 'psl'


if (grepl(var.name,'(pr|psl|huss)'))
  prct <- TRUE


gcm <- 'CanESM2'

##Load RCM Data
read.dir <- '/storage/data/climate/downscale/BCCAQ2/CanRCM4/'

rcm.dates <- seq(from=as.Date('1950-01-01'),by='day',to=as.Date('2100-12-31'))
flag <- grep('*-02-29',rcm.dates)
rcm.dates <- rcm.dates[-flag]

cells <- c(76,150)

rcm.series <- get.rcm.data(var.name,cells,read.dir)
rcm.seasonal <- get.rcm.seasonal(var.name,rcm.series,rcm.dates)

##---------------------------------------------------------------------
##Load GCM Data
cell.subset <- c(85,50) ##Nanaimo

data.dir <- paste('/storage/data/climate/downscale/BCCAQ2/CMIP5/',gcm,sep='')
gcm.data <- get.gcm.data(var.name,gcm,data.dir,cell.subset)

if (var.name=='huss') {
   rcm.series <- rcm.series*1000
   gcm.data$data <- gcm.data$data*1000
}


gcm.time.1980s <- get.data.subset(gcm.data$time,gcm.data$time,'1971-2000')
gcm.time.2050s <- get.data.subset(gcm.data$time,gcm.data$time,'2041-2070')

gcm.data.1980s <- get.data.subset(gcm.data$data,gcm.data$time,'1971-2000')
gcm.data.2020s <- get.data.subset(gcm.data$data,gcm.data$time,'2011-2040')
gcm.data.2050s <- get.data.subset(gcm.data$data,gcm.data$time,'2041-2070')
gcm.data.2080s <- get.data.subset(gcm.data$data,gcm.data$time,'2071-2100')

gcm.clim.1980s <- mean(gcm.data.1980s,na.rm=T)
gcm.clim.2050s <- mean(gcm.data.2050s,na.rm=T)
gcm.anoms <- gcm.clim.2050s - gcm.clim.1980s

gcm.seas.1980s <- get.season.data(gcm.data.1980s,gcm.time.1980s,'mean')
gcm.seas.2050s <- get.season.data(gcm.data.2050s,gcm.time.2050s,'mean')

gcm.seasmax.1980s <- get.season.data(gcm.data.1980s,gcm.time.1980s,'max')
gcm.seasmax.2050s <- get.season.data(gcm.data.2050s,gcm.time.2050s,'max')

gcm.seas.anoms <- gcm.seas.2050s - gcm.seas.1980s
gcm.seasmax.anoms <- (gcm.seasmax.2050s - gcm.seasmax.1980s)/gcm.seasmax.1980s*100
if (prct) {
  gcm.seas.anoms <- (gcm.seas.2050s - gcm.seas.1980s)/gcm.seas.1980s*100
  gcm.seasmax.anoms <- gcm.seasmax.2050s - gcm.seasmax.1980s
}

##**********************************************************************************
##**********************************************************************************
##Plotting
var.range <- list(psl=c(101.2,102.2),
                  pr=c(1000,3000),
                  huss=c(4,12),
                  tas=c(4,18),
                  uas=c(0,2),   
                  vas=c(0,2.25),
                  wspd=c(2,5))
var.seq <- list(psl=seq(101,102.2,0.2),
                pr=seq(1000,3000,500),
                huss=seq(4,12,2),
                tas=seq(4,18,2),
                uas=seq(0,2,0.25),      
                vas=seq(0,2.25,0.25),
                wspd=seq(2,5,0.5))
range.max <- list(psl=c(103.1,105.1),
                  pr=c(30,120),
                  huss=c(10,24),
                  tas=c(15,45),
                  uas=c(5,10),   
                  vas=c(5,13),
                  wspd=c(6,14))
seq.max <- list(psl=seq(103.1,105.1,0.2),
                pr=seq(30,120,30),
                huss=seq(8,24,4),
                tas=seq(15,45,5),
                uas=seq(5,10,1),      
                vas=seq(5,13,1),
                wspd=seq(6,14,2))

var.units <- list(psl='kPa',
                pr='mm',
                huss='kg/kg',
                tas='degC',
                uas='m/s',      
                vas='m/s',
                wspd='m/s')


   yr.rcm.fac <- as.factor(format(rcm.dates,'%Y'))
   yrs <- as.numeric(levels(yr.rcm.fac))
   rcm.ann.mean <- tapply(rcm.series,yr.rcm.fac,max)
   print(range(rcm.ann.mean))
   yr.fac <- as.factor(format(as.Date(gcm.data$time),'%Y'))
   ann.avg <- tapply(gcm.data$data,yr.fac,max,na.rm=T)
   print(range(ann.avg))
 

plot.dir <- '/storage/data/projects/rci/building_code/gcm_plots/'

compare.annual.series(gcm.data,rcm.series,rcm.dates,plot.dir,
                      var.name,var.units[[var.name]],fxn='mean',var.range[[var.name]],var.seq[[var.name]])
                      
compare.annual.series(gcm.data,rcm.series,rcm.dates,plot.dir,
                      var.name,var.units[[var.name]],fxn='max',range.max[[var.name]],seq.max[[var.name]])

##***********************************************************
##Windspeed
if (1==1) {
uas.rcm.series <- get.rcm.data('uas',cells,read.dir)
vas.rcm.series <- get.rcm.data('vas',cells,read.dir)
wspd.rcm.series <- sqrt(uas.rcm.series^2 + vas.rcm.series^2)

uas.gcm.data <- get.gcm.data('uas',gcm,data.dir,cell.subset)
vas.gcm.data <- get.gcm.data('vas',gcm,data.dir,cell.subset)
wspd.gcm.data <- uas.gcm.data
wspd.gcm.data$data <- sqrt(uas.gcm.data$data^2 + vas.gcm.data$data^2)

compare.annual.series(wspd.gcm.data,wspd.rcm.series,rcm.dates,plot.dir,
                      'wspd',var.units[['wspd']],fxn='mean',var.range[['wspd']],var.seq[['wspd']])                      

compare.annual.series(wspd.gcm.data,wspd.rcm.series,rcm.dates,plot.dir,
                      'wspd',var.units[['wspd']],fxn='max',range.max[['wspd']],seq.max[['wspd']])
}
