##Script to plot the fraction of useable MODIS data
library(raster)
library(rgdal)
library(TeachingDemos)
library(maps)
library(graticule)

library(ncdf4)

##X-Axis Ticks
get.proj.xaxis <- function(lons,crs,plot.window.ylim) {

  y <- seq(0,80,0.1)
  xm <- sapply(lons,rep,length(y))
  S <- apply(xm,2,function(x) {SpatialPoints(cbind(x,y), proj4string = CRS("+proj=longlat +datum=WGS84"))})
  S2<- lapply(S,spTransform, crs)
  indices <- lapply(S2,function(x){which.min(abs(x@coords[,'y']-plot.window.ylim[1]))})
  xticks <- mapply(FUN=function(s,i){s@coords[,'x'][i]},S2,indices)
  return(xticks)
}

 ##Y-Axis Ticks
get.proj.yaxis <- function(lats,crs,plot.window.xlim) {

  x <- seq(-180,-80,0.1)
  ym <- sapply(lats,rep,length(x))
  S <- apply(ym,2,function(y) {SpatialPoints(cbind(x,y), proj4string = CRS("+proj=longlat +datum=WGS84"))})
  S2<- lapply(S,spTransform, crs)
  indices <- lapply(S2,function(x){which.min(abs(x@coords[,'x']-plot.window.xlim[1]))})
  yticks <- mapply(FUN=function(s,i){s@coords[,'y'][i]},S2,indices)
  return(yticks)
}



make.plot.window <- function(bounds) {

  xleft  <- 0.1
  xright <- -0.05
  ybot   <- 0.05
  ytop   <- -0.05

  xlim.min <- bounds@xmin
  xlim.max <- bounds@xmax
  ylim.min <- bounds@ymin
  ylim.max <- bounds@ymax

  ##Set plot boundaries
  xlim.adj <- (xlim.max - xlim.min)
  ylim.adj <- (ylim.max - ylim.min)
  plot.window.xlim <- c((xlim.min + xlim.adj*xleft), (xlim.max + xlim.adj*xright))
  plot.window.ylim <- c((ylim.min + ylim.adj*ybot),  (ylim.max + ylim.adj*ytop))
  rv <- list(xlim=plot.window.xlim,
             ylim=plot.window.ylim)
  return(rv)
}

get.region.shape <- function(region,shape.dir) {
  region.shp <- readOGR(shape.dir, region, stringsAsFactors=F, verbose=F)
  return(region.shp)
}

add.graticules <- function(lons,lats,crs) {

  xl <-  range(lons)
  yl <- range(lats) 
  ###grat <- graticule(lons, lats, proj = CRS(crs),xlim=xl,ylim=xl)
  grat <- graticule(lons, lats, proj = CRS(crs),xlim=c(-145,-110),ylim=c(48,62))
  ###labs <- graticule_labels(lons = lons, lats = lats, xline=lons,yline=lats,proj = CRS(crs))
  ##labs <- graticule_labels(lons = lons, lats = lats, xline=lons[2],yline=lats[2],proj = CRS(crs))
  
  ##rv <- list(grat=grat,labs=labs,lons=lons,lats=lats)
  rv <- grat
  return(rv)
}
        

alb.crs <- "+init=epsg:3005"
lons <- c(-145.0, -140.0,-135.0,-130.0,-125.0,-120.0,-115.0,-110.0)
lats <- c(  48.0,  50.0,  52.0,  54.0,  56.0,  58.0,  60.0,  62.0)

grats <- add.graticules(lons,lats,alb.crs)

file.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/ERA5/concat/tps_era5/' 
file <- "pr.mon.sum.ERA5-tps.nc"
brick.data <- brick(paste0(file.dir,file))

mn <- 4

box.data <- subset(brick.data,mn)
box.data <-projectRaster(box.data,crs=CRS(alb.crs))
bounds <- extent(box.data)
plot.window <- make.plot.window(bounds)
box.range <- round(c(box.data@data@min,box.data@data@max)/100)*100
 
nc <- nc_open(paste0(file.dir,file))
lon <- ncvar_get(nc,'lon')
lat <- ncvar_get(nc,'lat')
nc_close(nc)

colour.map <- function() {
  blue.to.light.blue <- colorRampPalette(colors=c("#0000CC","#33FFFF"),bias=1, space = "Lab", interpolate = "linear")
  light.blue.to.green <- colorRampPalette(colors=c("#33FFFF",'#99FF99'),bias=1, space = "Lab", interpolate = "linear")
  green.to.tan <- colorRampPalette(colors=c('#99FF99','#FFDD97'),bias=1, space = "Lab", interpolate = "linear")
  tan.to.brown <- colorRampPalette(colors=c('#FFDD97','#FFFF99','#C1934D'),bias=1, space = "Lab", interpolate = "linear")

  rv <- c(blue.to.light.blue(11),light.blue.to.green(20),green.to.tan(20),tan.to.brown(11))
  return(rv)
}


##---------------------------------------------------------------------------------
##Make Plot

  region <- 'bc'
  plot.file <- paste0('/storage/data/projects/rci/building_code/ERA5/',
                      var.name,'.',month.abb[mn],'.ERA5-TPS.png')
  shape.dir <- paste('/storage/data/projects/rci/data/assessments/shapefiles/bc/',sep='')
  region.shp <- spTransform(get.region.shape(region,shape.dir),CRS(alb.crs))

  ##width <- 1000 ##3926
  ##height <- 900 ##2383
  png(file=plot.file,width=8,height=5,units='in',res=600,pointsize=6,bg='white')
  layout(matrix(c(1,1,1,1,1,1,1,1,2,3),nrow=2))
  par(mar=c(6,6,6,1))
  plot(c(),xlim=plot.window$xlim,ylim=plot.window$ylim,xaxs='i',yaxs='i',
     bg='lightgray',axes=FALSE,
     xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main='',
     cex.axis=2.5,cex.lab=2.5,cex.main=2.4,mgp=c(3.5,2,0))
     rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='white')
     glen <- dim(grats$labs@coords)[1]
     ghalf <- glen/2

     xtks <- get.proj.xaxis(lons,alb.crs,plot.window$ylim)
     ytks <- get.proj.yaxis(lats,alb.crs,plot.window$xlim)
     axis(2,at=ytks,label=lats,cex.axis=1.75)
     axis(1,at=xtks,label=lons,cex.axis=1.75)

     ###axis(1,at=unclass(grats$labs@coords)[1:ghalf,1],label=grats$lons,cex.axis=1.75)
     ###axis(2,at=unclass(grats$labs@coords)[(ghalf+1):glen,2],label=grats$lats,cex.axis=1.75)

     bc.dir <- '/storage/data/gis/basedata/base_layers/'
     bc.overlay <- 'bc_province_wgs84'
     bc.shp <- readOGR(bc.dir, bc.overlay, stringsAsFactors=F, verbose=F)
     us.shp <- readOGR(bc.dir, 'united_states', stringsAsFactors=F, verbose=F)

     image(box.data,add=T,col = colour.map())    

     plot(spTransform(us.shp,CRS(alb.crs)),add=TRUE,border='black',lwd=1)
     plot(spTransform(bc.shp,CRS(alb.crs)),add=TRUE,border='black',lwd=1)
     plot(grats$grat,add=TRUE,lty=3,col='gray',lwd=1)

     box(which='plot',lwd=2)
 
u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)     

##    legend('bottomright',legend=c(col=colour.map(),cex=2,pch=15,title='Pr Diff (mm)')


v <- c( 0.8*v[2], v[2], 0.75*v[4], v[4] )
##v <- c( v[2], 1.2*v[2], 0.75*v[4], v[4] )
     bc.proj <- spTransform(bc.shp,CRS("+init=epsg:3005"))  
     pnw.bnds <- spTransform(readOGR('/storage/data/projects/rci/data/winter_sports/study_map','north_america_state_provincial_boundaries', stringsAsFactors=F, verbose=F),CRS("+init=epsg:3005"))  
     bounds <- extent(bc.proj)
     xlim <- c(bounds@xmin,bounds@xmax)
     ylim <- c(bounds@ymin,bounds@ymax)

##par( fig=v, new=TRUE, mar=c(0,0,0,0) )
  par(mar=c(10,0,6,3))
  plot(c(),xlim=xlim,ylim=ylim,xaxs='i',yaxs='i',
     bg='gray94',# 'lightgray',
     xlab='',ylab='',main='',axes=F)
     rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='lightblue')
     plot(pnw.bnds,add=T,col='darkgray')
     text(0.6*xlim[2],0.6*ylim[2],'British\nColumbia',cex=1.7)
     box(which='plot',lwd=2)

legend_image <- as.raster(matrix(colour.map(), ncol=1))
  par(mar=c(1,0,1,1))
  plot(c(0,3),c(0,3),type = 'n', axes = F,xlab = '', ylab = '', main='',yaxs='i',xaxs='i')
  rect(0.01,0.35,1,2.76,border='black',lwd=2)
  text(x=0.5, y = 2.9, labels ='Elev. (m)',cex=2.75)
  text(x=1.3, y = c(0.4,1.5,2.7), labels =c(box.range[1],(box.range[2]+box.range[1])/2,box.range[2]),cex=2.5)
  rasterImage(legend_image, 0.03, 0.36, 0.98,2.75)
##  box(which='plot')
  dev.off()