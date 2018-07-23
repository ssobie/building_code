#Create windspeed from uas and vas files
##'ACCESS1-0',
gcm.list <- c('CanESM2',
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
##gcm.list <- 'ACCESS1-0'
base.dir <- '/storage/data/climate/downscale/CMIP5/building_code/'

for (gcm in gcm.list) {
  ##Windspeed
  read.dir <- paste0(base.dir,gcm,'/')
  write.dir <- read.dir
  all.files <- list.files(path=read.dir,pattern='rcp85')
  uas.file <- all.files[grep('uas',all.files)]
  vas.file <- all.files[grep('uas',all.files)]
  wspd.out <- gsub('uas','wspd',uas.file)

  uas2 <- paste0('cdo -O mul ',read.dir,uas.file,' ',read.dir,uas.file,' ',write.dir,'uas2.nc')
  print(uas2)
  system(uas2)
  vas2 <- paste0('cdo -O mul ',read.dir,vas.file,' ',read.dir,vas.file,' ',write.dir,'vas2.nc')
  print(vas2)
  system(vas2)

  wspd2 <- paste0('cdo -O add ',write.dir,'uas2.nc ',write.dir,'vas2.nc ',write.dir,'wspd2.nc')
  print(wspd2)
  system(wspd2)

  ##Windspeed
  wspd <- paste0('cdo -O sqrt ',write.dir,'wspd2.nc ',read.dir,wspd.out)
  print(wspd)
  system(wspd)

  work <- paste0('ncrename -v uas,wspd ',read.dir,wspd.out)
  print(work)
  system(work)

  file.remove(paste0(read.dir,'uas2.nc'))
  file.remove(paste0(read.dir,'vas2.nc'))
  file.remove(paste0(read.dir,'wspd2.nc'))
}