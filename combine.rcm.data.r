##Script to combine the separate year files from CanRCM4,
##stiching together individual grid cells

library(ncdf4)
library(rgdal)
library(raster)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

rcm.dir <- '/storage/data/climate/downscale/RCM/CanRCM4/daily/NAM_22/'
write.dir <- '/storage/data/climate/downscale/RCM/CanRCM4/daily/'
var.list <- c('pr','tas','tasmax','tasmin','uas','vas','ps','huss','snd','snc')
var.list <- 'snd'


for (var.name in var.list) {
print(var.name)
##Nanaimo
##cells <- c(76,150)
##lon.i <- -123.968641 
##lat.i <- 49.185618

##Cowichan Hospital
##lon.i <- -123.722256
##lat.i <- 48.785912

##Granville and 41st Hospital
lon.i <- -123.139690
lat.i <- 49.234305

##---------------------------------------------
##Find the closest grid cell to a point
initial.file <- paste(rcm.dir,var.name,'/',var.name,'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-19501231.nc',sep='')
hist.prefix <- paste(var.name,'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_',sep='')
proj.prefix <- paste(var.name,'_NAM-22_CCCma-CanESM2_rcp85_r1i1p1_CCCma-CanRCM4_r2_day_',sep='')
data <- raster(initial.file)

nc <- nc_open(initial.file)
lon <- ncvar_get(nc,'lon')
lon <- ((lon + 180) %% 360) - 180
lat <- ncvar_get(nc,'lat')

regular <- "+proj=longlat +ellps=WGS84"
rotated <- "+proj=ob_tran +o_proj=longlat +lon_0=-97 +o_lat_p=42.5 +a=1 +to_meter=0.0174532925199 +no_defs"

projection(data) <- rotated
my.extent <- extent(data)
xy.rot <- coordinates(data)
sp.rot <- SpatialPoints(xy.rot)
projection(sp.rot) <- rotated
sp.reg <- spTransform(sp.rot,CRS(regular))
xy.reg <- coordinates(sp.reg)

dat <- as.data.frame(xy.reg[,2])
colnames(dat) <- "lat"
spdf <- SpatialPointsDataFrame(sp.rot,dat)
dat <- rasterize(spdf,data)
dat <- subset(dat,2:2)

geod <- sqrt((lon - lon.i)^2 + (lat - lat.i)^2)
coord.ix <- which(geod==min(geod),arr.ind=TRUE)
##cells <- as.numeric(coord.ix)
cells <- c(80,152)
coord.ix <- cells

##----------------------------------------------
hist.ysts <- seq(1951,2001,5)
hist.yens <- seq(1955,2005,5)

proj.ysts <- seq(2006,2096,5)
proj.yens <- seq(2010,2100,5)

past.years <- paste(hist.ysts,'0101-',hist.yens,'1231.nc',sep='')
past.list <- paste(rcm.dir,var.name,'/',hist.prefix,past.years,sep='')

proj.years <- paste(proj.ysts,'0101-',proj.yens,'1231.nc',sep='')
proj.list <- paste(rcm.dir,var.name,'/',proj.prefix,proj.years,sep='')

file.list <- c(initial.file,past.list,proj.list)

##Loop over the files, grab the grid cell and stitch
##them together as one file

flen <- length(file.list)
data <- c()
time <- c()

for (i in seq_along(file.list)) {
    file <- file.list[i]
    print(file)
    fnc <- nc_open(file)
    data.sub <- ncvar_get(fnc,start=c(coord.ix[1],coord.ix[2],1),count=c(1,1,-1))
    data <- c(data,data.sub)
    sub.time <- as.character(format(netcdf.calendar(fnc),'%Y-%m-%d'))     
    time <- c(time,sub.time)

    nc_close(fnc)        
}

if (var.name=='pr') {
   data.final <- data*86400
   data.final[data.final<0] <- 0
   data <- data.final
}
if (grepl('(tasmax|tasmin)',var.name)) {
   data.final <- data - 273
   data <- data.final
}


rcm.data <- list(data=data,time=time)

save(rcm.data,file=paste(write.dir,var.name,'_GRANVILLE_41st_',cells[1],'_',cells[2],'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))


}

