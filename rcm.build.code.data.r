##Contains functions to compute Building Code values from the BC
##Building Code and the American ASHRAE Handbook


##************************************************************************
##************************************************************************
##Load Data
read.dir <- '/storage/data/climate/downscale/BCCAQ2/CanRCM4/'

dates <- seq(from=as.Date('1950-01-01'),by='day',to=as.Date('2100-12-31'))
flag <- grep('*-02-29',dates)
dates <- dates[-flag]
cell <- '76_150'
load(paste(read.dir,'tas_NANAIMO_',cell,'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
tas.data <- data - 273

load(paste(read.dir,'tasmax_NANAIMO_',cell,'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
tasmax.data <- data

load(paste(read.dir,'tasmin_NANAIMO_',cell,'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
tasmin.data <- data

load(paste(read.dir,'pr_NANAIMO_',cell,'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
pr.data <- data
pr.data[pr.data < 0.1] <- 0

load(paste(read.dir,'ps_NANAIMO_',cell,'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
pas.data <- data/1000

load(paste(read.dir,'huss_NANAIMO_',cell,'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
huss.data <- data

load(paste(read.dir,'uas_NANAIMO_',cell,'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
uas.data <- data

load(paste(read.dir,'vas_NANAIMO_',cell,'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
vas.data <- data

load(paste(read.dir,'snd_NANAIMO_',cell,'_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_19500101-21001231.RData',sep=''))
snd.data <- data


