##Script to combine the separate year files from CanRCM4,
##stiching together individual grid cells

library(ncdf4)
##library(rgdal)
library(raster)
library(PCICt)

##  args <- commandArgs(trailingOnly=TRUE)
##  for(i in 1:length(args)){
##      eval(parse(text=args[[i]]))
##  }
var.name <- 'vas'
print(var.name)
rcm.dir <- paste0('/storage/data/climate/downscale/RCM/CanRCM4/1hr/atmos/run1/',var.name,'/r10i2p1/')
write.dir <- '/storage/data/climate/downscale/BCCAQ2/CanRCM4/'

cells <- c(39,77)

lon.i <- -123.968641 
lat.i <- 49.185618

##---------------------------------------------
##Find the closest grid cell to a point
initial.file <- list.files(rcm.dir,pattern='*.nc',full.name=TRUE)[1]

##raw <- raster(initial.file)
##data <- subset(raw,1)
##rm(raw)

nc <- nc_open(initial.file)
lon <- ncvar_get(nc,'lon')
lon <- ((lon + 180) %% 360) - 180
lat <- ncvar_get(nc,'lat')

if (1==0) {
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
}

coord.ix <- cells

##----------------------------------------------
ysts <- seq(1950,2099,1)
yens <- seq(1951,2100,1)

file.list <- list.files(rcm.dir,pattern='*.nc',full.name=TRUE)


##Loop over the files, grab the grid cell and stitch
##them together as one file

flen <- length(file.list)
data <- c()
series <- c()

for (i in seq_along(file.list)) {
    file <- file.list[i]
    print(file)
    fnc <- nc_open(file)
    data.sub <- ncvar_get(fnc,start=c(coord.ix[1],coord.ix[2],1),count=c(1,1,-1))
    time.atts <- ncatt_get(fnc,'time')
    time.calendar <- time.atts$calendar
    time.units <- time.atts$units
    time.values <- ncvar_get(fnc,'time')
    origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                             cal=time.calendar)
    time.series <- (origin.pcict + time.values*86400)+1
    time.chars <- as.character(time.series)
    data <- c(data,data.sub)
    series <- c(series,time.chars)
    print(range(time.series))
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

series <- gsub(':00:01',':00:00',series)

output <- list(data=data,time=series)

save(output,file=paste0("/storage/data/climate/downscale/RCM/CanRCM4/van/",var.name,'_vancouver_NAM-44_CCCma-CanESM2_historical-r2_r10i2p1_CCCma-CanRCM4_r2_1hr_19500101-21001231.nc'))




