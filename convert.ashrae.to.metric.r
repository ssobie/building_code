##Script to convert the ASHRAE Parameters from Imperial to SI Units

celsius <- function(x){(x-32)/1.8} ##F to C
pr.mm <- function(x){x*25} ##inches to mm
enthap.si <- function(x){x*2.326} ##Btu/lb to kJ/kg
wspd.si <- function(x){x*1.6093/3.6} ##mph to m/s
hr.si <- function(x){x/7} ##grains/lb to g/kg
deg.si <- function(x){x} ##preserve degrees

convert.to.si <- function(x) {
  data <- as.numeric(x[1])
  unit <- x[2]  
  unit.fxn <- switch(unit,
                     degF=celsius,
                     ratio=hr.si,
                     btu.lb=enthap.si,
                     mph=wspd.si,
                     inch=pr.mm,        
                     deg=deg.si)   
  si <- unit.fxn(data)                                     
  return(si)
}

##Nanaimo Variables
ashrae.vars <-  list(c(11.2,'degF'), ##Cold Dewpoint Temp 0.4% degF
                     c(9.8,'ratio'),  ##Cold Humidity Ratio 0.4% gr/lb
                     c(28.4,'degF'), ##MCDB Dewpoint 0.4% degF
                     c(63.9,'degF'), ##MCWB Dry Bulb 99.6% degF
                     c(65.2,'degF'), ##Evap Wet Bulb 99.6% degF
                     c(76.8,'degF'), ##MCDB Wet Bulb 99.6% degF
                     c(60.8,'degF'), ##Warm Dewpoint 99.6% degF (Dehumidification)
                     c(79.9,'ratio'), ##Warm Humidity Ratio 99.6% gr/lb 
                     c(68.9,'degF'),  ##MCDB Depoint 99.6% degF
                     c(30.1,'btu.lb'), ##Enthalpy 99.6% Btu/lb
                     c(77.0,'degF'),  ##MCDB Enthalpy 99.6% degF
                     c(75.2,'degF'), ##Annual Max Wet Bulb degF
                     c(23.4,'degF'), ##Heating Dry Bulb 0.4% degF
                     c(80.3,'degF'), ##Cooling Dry Bulb 99.6% degF
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
                     c(25.9,'mph'),  ##Extreme Annual Wind Speed 97.5% mph
                     c(30.2,'mph'),  ##Extreme Annual Wind Speed 99% mph
                     c(89.9,'degF'), ##5-Year Return Period Max Temp degF
                     c(17.0,'degF'), ##5-Year Return Period Min Temp degF
                     c(NA,'inch')) ##5-Year Return Period Daily Precip inches
                     
