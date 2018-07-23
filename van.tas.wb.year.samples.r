##Script to create the year-long series of hourly temperature and wet bulb temp for Vancouver
##2000,2020,2050,2080

source('/storage/home/ssobie/code/repos/building_code/building.code.fxns.r',chdir=TRUE)

read.dir <- '/storage/data/climate/downscale/RCM/CanRCM4/van/'

load(paste0(read.dir,'huss_vancouver_NAM-44_CCCma-CanESM2_historical-r2_r10i2p1_CCCma-CanRCM4_r2_1hr_19500101-21001231.nc'))
huss.output <- output
rm(output)

load(paste0(read.dir,'psl_vancouver_NAM-44_CCCma-CanESM2_historical-r2_r10i2p1_CCCma-CanRCM4_r2_1hr_19500101-21001231.nc'))
psl.output <- output
rm(output)

load(paste0(read.dir,'tas_vancouver_NAM-44_CCCma-CanESM2_historical-r2_r10i2p1_CCCma-CanRCM4_r2_1hr_19500101-21001231.nc'))
tas.output <- output
rm(output)

##Dew Point Temperature
dpt <- dew.point.temp(psl.output$data/1000,huss.output$data) 

##Wet Bulb Temperature
wbt <- temp.wet.bulb(tas.output$data,dpt,psl.output$data,huss.output$data)

intervals <- c('2000','2020','2050','2080')

time.list <- vector(mode='list',length=4)
tas.list <- vector(mode='list',length=4)
wbt.list <- vector(mode='list',length=4)

for (i in seq_along(intervals)) {
  ix <- grep(paste0(intervals[i],'-*'),tas.output$time)
  time.list[[i]] <- tas.output$time[ix]
  tas.list[[i]] <- tas.output$data[ix]
  wbt.list[[i]] <- wbt[ix]
}

par(mfrow=c(2,1))
plot(c(),xlim=c(0,8770),ylim=range(unlist(tas.list)))
lines(tas.list[[1]],col='black')
lines(tas.list[[2]],col='yellow')
lines(tas.list[[3]],col='orange')
lines(tas.list[[4]],col='red')

plot(c(),xlim=c(0,8770),ylim=range(unlist(wbt.list)))
lines(wbt.list[[1]],col='black')
lines(wbt.list[[2]],col='yellow')
lines(wbt.list[[3]],col='orange')
lines(wbt.list[[4]],col='red')


data.out <- rbind(c('Time','TAS','WetBulb'),
                  cbind(unlist(time.list),round(unlist(tas.list)-273,1),round(unlist(wbt.list)-273)))

write.table(data.out,file=paste0(read.dir,'van_tas_wbt_year_subsets.csv'),row.name=F,col.name=F,quote=F,sep=',')