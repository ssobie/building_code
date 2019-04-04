##Script to convert the downloaded ERA5 data into useful
##file formats and structures

library(ncdf4)



read.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/ERA5/concat/sim/'
tmp.dir <- '/local_temp/ssobie/era5/'
if (!file.exists(tmp.dir))
   dir.create(tmp.dir,recursive=TRUE)

##-------------------------------------------------------------------------------------
##Make tas - tasmin
file.list <- list(tas='tas_QDM_day_ERA5_BC_19800101-20181231.nc',
                  tasmax='tasmax_QDM_day_ERA5_BC_19800101-20181231.nc',
                  tasmin='tasmin_QDM_day_ERA5_BC_19800101-20181231.nc',
                  tasrange='tasrange_QDM_day_ERA5_BC_19800101-20181231.nc',
                  tasskew='tasskew_QDM_day_ERA5_BC_19800101-20181231.nc')

tasmax.alt.file <- 'alt_tasmax_QDM_day_ERA5_BC_19800101-20181231.nc'
tasmin.alt.file <- 'alt_tasmin_QDM_day_ERA5_BC_19800101-20181231.nc'

for (file in file.list) {
  file.copy(from=paste0(read.dir,file),to=tmp.dir,overwrite=TRUE)
}
file.copy(from=paste0(tmp.dir,file.list$tasmax),to=paste0(tmp.dir,tasmax.alt.file),overwrite=TRUE)
file.copy(from=paste0(tmp.dir,file.list$tasmin),to=paste0(tmp.dir,tasmin.alt.file),overwrite=TRUE)


tas.nc <- nc_open(paste0(tmp.dir,file.list$tas))
tasrange.nc <- nc_open(paste0(tmp.dir,file.list$tasrange)) 
tasskew.nc <- nc_open(paste0(tmp.dir,file.list$tasskew))
alttasmax.nc <- nc_open(paste0(tmp.dir,tasmax.alt.file),write=TRUE)
alttasmin.nc <- nc_open(paste0(tmp.dir,tasmin.alt.file),write=TRUE)

tas <- ncvar_get(tas.nc,'tas')
tasskew <- ncvar_get(tasskew.nc,'tasskew')
tasrange <- ncvar_get(tasrange.nc,'tasrange')

tas.diff <- tasskew * tasrange
alttasmin <- tas - tas.diff
alttasmax <- alttasmin + tasrange

ncvar_put(alttasmax.nc,'tasmax',alttasmax)
ncvar_put(alttasmin.nc,'tasmin',alttasmin)

nc_close(alttasmax.nc)
nc_close(alttasmin.nc)

file.copy(from=paste0(tmp.dir,tasmax.alt.file),to=read.dir,overwrite=TRUE)
file.copy(from=paste0(tmp.dir,tasmin.alt.file),to=read.dir,overwrite=TRUE)

nc_close(tas.nc)
nc_close(tasrange.nc)
nc_close(tasskew.nc)

for (file in file.list) {
  file.remove(paste0(tmp.dir,file))
}
file.remove(paste0(tmp.dir,tasmax.alt.file))
file.remove(paste0(tmp.dir,tasmin.alt.file))
