##Script to contain the comments for the openxlsx summary tables


variable.comment <- function(var.name) {

   var.list <- list(pr.ann.avg=c('Precipitation\n','Average annual precipitation total'),
                    pr.ann.max=c('Max Ann Pr\n','Maximum annual precipitation amount'),
                    pr_rp5=c('RP5 Precipitation\n','5-Year annual maximum one day precipitation amount'),
                    pr_rp20=c('RP20 Precipitation\n','20-Year annual maximum one day precipitation amount'),
                    pr_rp50=c('RP50 Precipitation\n','50-Year annual maximum one day precipitation amount'),
                    cdd.18.ann=c('CDD 18.0 degC\n','Cooling Degree Days (Threshold: 18.0\u00B0C)'),
                    cdd.18.3.ann=c('CDD 18.3 degC\n','Cooling Degree Days (Threshold: 18.3\u00B0C)'),
                    hdd.18.ann=c('HDD 18.0 degC\n','Heating Degree Days (Threshold: 18.0\u00B0C)'),
                    hdd.18.3.ann=c('HDD 18.3 degC\n','Heating Degree Days (Threshold: 18.3\u00B0C)'),
                    tasmax_rp5=c('RP5 Tasmax\n','5-Year annual maximum daily maximum temperature'),
                    tasmax_rp20=c('RP20 Tasmax\n','20-Year annual maximum daily maximum temperature'),
                    tasmin_rp5=c('RP5 Tasmin\n','20-Year annual minimum daily minimum temperature'),
                    tasmin_rp20=c('RP20 Tasmin\n','20-Year annual minimum daily minimum temperature'),
                    tasmax.950.mon=c('Tasmax 95.0%\n','Warm Month Design Temperature 95.0%'),
                    tasmax.975.mon=c('Tasmax 97.5%\n','Warm Month Design Temperature 97.5%'),
                    tasmax.990.mon=c('Tasmax 99.0%\n','Warm Month Design Temperature 99.0%'),
                    tasmax.996.mon=c('Tasmax 99.6%\n','Warm Month Design Temperature 99.6%'),
                    tasmin.004.mon=c('Tasmin 0.4%\n','Cold Month Design Temperature 0.4%'),
                    tasmin.010.mon=c('Tasmin 1.0%\n','Cold Month Design Temperature 1.0%'),
                    tasmin.025.mon=c('Tasmin 2.5%\n','Cold Month Design Temperature 2.5%'),
                    tasmin.050.mon=c('Tasmin 5.0%\n','Cold Month Design Temperature 5.0%'),
                    wb.950.mon=c('Wetbulb 95.0%\n','Warm Month Wetbulb Temperature 95.0%'),
                    wb.975.mon=c('Wetbulb 97.5%\n','Warm Month Wetbulb Temperature 97.5%'),
                    wb.990.mon=c('Wetbulb 99.0%\n','Warm Month Wetbulb Temperature 99.0%'),
                    wb.996.mon=c('Wetbulb 99.6%\n','Warm Month Wetbulb Temperature 99.6%'),
                    snow_rp5=c('RP5 Snow\n','5-Year annual maximum daily snow load'),
                    snow_rp50=c('RP50 Snow\n','50-Year annual maximum daily snow load'),
                    wind_rp5=c('RP5 Wind\n','5-Year annual maximum daily windspeed'),
                    wind_rp10=c('RP10 Wind\n','10-Year annual maximum daily windspeed'),
                    wind_rp50=c('RP50 Wind\n','50-Year annual maximum daily windspeed'))


  rv <- var.list[[var.name]]
  return(rv)                    
}