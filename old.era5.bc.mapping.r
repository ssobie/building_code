##Script to map ERA5 data

library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(ncdf4)

source('/storage/home/ssobie/code/repos/assessments/resource.region.map.support.r',chdir=T)
source('/storage/home/ssobie/code/repos/assessments/bc_albers_map_support.r',chdir=T)

source(paste0('/storage/home/ssobie/code/repos/assessments/bc_map_support.r'),chdir=T)

get.region.shape <- function(region,shape.dir) {
  region.shp <- readOGR(shape.dir, region, stringsAsFactors=F, verbose=F)
  return(region.shp)
}


region <- 'bc'
meta <- get.region.names(region)
shape.dir <- paste('/storage/data/projects/rci/data/assessments/shapefiles/',meta$subset,'/',sep='')
region.shp <- spTransform(get.region.shape(region,shape.dir),CRS("+init=epsg:4326")) ##Keep this projection to extract the data from lat/lon
leg.loc <- get.leg.loc(region)

var.name <- 'pr'

file <- "/storage/data/climate/downscale/BCCAQ2+PRISM/ERA5/concat/tps_era5/pr.mon.sum.ERA5-tps.nc"

box <- brick(file)
for (mn in 1:12) {
  print(mn)
  mon.box <- subset(box,mn)
  region.range <- box.range <- range(as.matrix(mon.box),na.rm=T)

  past.plot.file <- paste0('/storage/data/projects/rci/building_code/ERA5/',var.name,'.',month.abb[mn],'.ERA5-TPS.png')
  past.plot.title <- paste0('ERA5-TPS ',month.abb[mn],' Total Precipitation Difference (1980-2012)')

  reg.ds.maps(mon.box,region,region.range,box.range,
            var.name,type='past','ERA5',region.shp,
              past.plot.file,past.plot.title,
              add.overlays=add.overlays,
              add.cities=add.cities,add.districts=add.districts,
              add.graticules=add.graticules,leg.loc=leg.loc,draft=FALSE)
}