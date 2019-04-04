##Building Code and the American ASHRAE Handbook

##------------------------------------------------------------------
##Nanaimo Variables - Don't have access to ASHARE 
ashrae.vars <-  list(c(11.2,'degF'), ##Cold Dewpoint Temp 0.4% degF
                     c(9.8,'ratio'),  ##Cold Humidity Ratio 0.4% gr/lb
                     c(28.4,'degF'), ##MCDB Dewpoint 0.4% degF
                     c(63.9,'degF'), ##MCWB Dry Bulb 99.6% degF
                     c(62.6,'degF'), ##MCWB Dry Bulb 99.0% degF
                     c(61.4,'degF'), ##MCWB Dry Bulb 98.0% degF
                     c(65.2,'degF'), ##Evap Wet Bulb 99.6% degF
                     c(63.8,'degF'), ##Evap Wet Bulb 99.0% degF
                     c(76.8,'degF'), ##MCDB Wet Bulb 99.6% degF
                     c(73.9,'degF'), ##MCDB Wet Bulb 99.0% degF
                     c(60.8,'degF'), ##Warm Dewpoint 99.6% degF (Dehumidification)
                     c(59.6,'degF'), ##Warm Dewpoint 99.0% degF (Dehumidification)
                     c(79.9,'ratio'), ##Warm Humidity Ratio 99.6% gr/lb
                     c(76.3,'ratio'), ##Warm Humidity Ratio 99.0% gr/lb
                     c(68.9,'degF'),  ##MCDB Depoint 99.6% degF
                     c(67.4,'degF'),  ##MCDB Depoint 99.0% degF
                     c(30.1,'btu.lb'), ##Enthalpy 99.6% Btu/lb
                     c(77.0,'degF'),  ##MCDB Enthalpy 99.6% degF
                     c(75.2,'degF'), ##Annual Max Wet Bulb degF
                     c(23.4,'degF'), ##Heating Dry Bulb 0.4% degF
                     c(27.1,'degF'), ##Heating Dry Bulb 1.0% degF                     
                     c(80.3,'degF'), ##Cooling Dry Bulb 99.6% degF
                     c(76.5,'degF'), ##Cooling Dry Bulb 99.0% degF
                     c(87.5,'degF'), ##Annual Max Dry Bulb degF
                     c(20.1,'degF'), ##Annual Min Dry Bulb degF
                     c(46.4,'inch'), ##Annual Total Precip inches
                     c(62.4,'inch'), ##Annual Max Total Precip inches
                     c(28.1,'inch'), ##Annual Min Total Precip inches
                     c(8.0,'inch'),  ##Annual SD Total Precip inches
                     c(40.4,'mph'),  ##High Annual Wind Speed 99.6% mph
                     c(41.8,'degF'),  ##Windy MCDB Windy Speed 99.6% degF
                     c(6.8,'mph'),   ##Cold MCWS Dry Bulb 0.4% mph
                     c(300,'deg'),   ##Cold PCWD Dry Bulb 0.4% degrees from north
                     c(7.5,'mph'),   ##Warm MCWS Dry Bulb 99.6% mph
                     c(340,'deg'),   ##Warm PCWD Dry Bulb 99.6% degrees from north
                     c(21.9,'mph'),  ##Extreme Annual Wind Speed 95% mph
                     c(25.9,'mph'),  ##Extreme Annual Wind Speed 97.5% mph
                     c(30.2,'mph'),  ##Extreme Annual Wind Speed 99% mph
                     c(89.9,'degF'), ##5-Year Return Period Max Temp degF
                     c(17.0,'degF'), ##5-Year Return Period Min Temp degF
                     c(NA,'inch')) ##5-Year Return Period Daily Precip inches

ashrae.si.versions <- unlist(lapply(ashrae.vars,convert.to.si))

##Vancouver (Granville & 41st) Variables
bc.si.versions <-  c(NA, ##'degC'), ##Cold Month Design Temperature 5.0% degC,
                     -6.0, ##'degC'), ##Cold Month Design Temperature 2.5% degC
                     -8.0, ##'degC'),  ##Cold Month Design Temperature 1.0% degC
                     NA, ##'degC'),  ##Cold Month Design Temperature 0.4% degC
                     NA,   ##'degC'), ##Warm Month Design Temperature 99.6% degC
                     NA,   ##'degC'), ##Warm Month Design Temperature 99.0% degC
                     28,   ##'degC'), ##Warm Month Design Temperature 97.5% degC
                     NA,   ##'degC'), ##Warm Month Design Temperature 95.0% degC
                     NA,   ##'degC'), ##Warm Month Wet Bulb Temperature 99.6% degC
                     NA,   ##'degC'), ##Warm Month Wet Bulb Temperature 99.0% degC
                     20,   ##'degC'), ##Warm Month Wet Bulb Temperature 97.5% degC
                     NA,   ##'degC'), ##Warm Month Wet Bulb Temperature 95.0% degC
                     2925, ##'DD'), ##Degree Days below 18.0 degC
                     NA,   ##'DD'), ##Degree Days below 18.3 degC
                     NA,   ##'DD'), ##Degree Days above 18.3 degC
                     NA,   ##'DD'), ##Degree Days above 18.0 degC
                     NA,   ##'mm'), ##Annual Max Daily Precipitation mm
                     1325, ##'mm'), ##Annual Average Total Precipitation mm
                     NA,   ##'mm'), ##5-Year return period Daily Precipitation
                     NA,   ##'mm'), ##10-Year return period Daily Precipitation
                     NA,   ##'mm'), ##20-Year return period Daily Precipitation
                     107,   ##'mm'), ##50-Year return period Daily Precipitation
                     NA,  ##'Pa'),  ##5-Year return period Wind Load
                     350,  ##'Pa'))  ##10-Year return period Wind Load
                     450,  ##'Pa'),  ##50-Year return period Wind Load
                     NA,   ##'kPa'), ##20-Year Return Period Snow Load kPa
                     1.9)  ##,'kPa')) ##50-Year Return Period Snow Load kPa

##------------------------------------------------------------------
