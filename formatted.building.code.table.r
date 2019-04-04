##Script to sort out the data requested for Okanagan

source('/storage/home/ssobie/code/repos/building_code/summary.BCBC.table.comments.r',chdir=T)

library(openxlsx)

##------------------------------------------------------------- 
##------------------------------------------------------------- 

get.units.pane <- function(var.name) {
  leg.label <- NA
  if (grepl("(tas|wb)", var.name))
    leg.label <- c(rep('degC',7),rep('%',6))
  if (grepl("(pr|rx|r10|r20|r9|RP|sdii|prcptot)", var.name))
    leg.label <- c(rep('mm',7),rep('%',6))
  if (grepl("(dd)", var.name))
    leg.label <- c(rep('degree days',7),rep('%',6))
  if (grepl("(fdE|cddE|cdd90|cddmax|cwd|su|gsl|id|trE|su30|r95daysE|r99daysE)", var.name))
    leg.label <- c(rep('days',7),rep('%',6))
  if (grepl("(wind)", var.name))
    leg.label <- c(rep('Pa',7),rep('%',6))
  if (grepl("(snow)", var.name))
    leg.label <- c(rep('kPa',7),rep('%',6))

  return(leg.label)
} 

##------------------------------------------------------------- 
get.round.val <- function(var.name) {
  rd <- 0
  if (grepl("(dd)", var.name))
    rd <- 0    
  if (grepl("(tas|cdd|hdd|wind|snow|wb)", var.name))
    rd <- 1
  if (grepl("(pr|rx|r9|RP|rp|tx90|tn10)", var.name))
    rd <- 0
  if (grepl("(pas|snowdepth|snow)", var.name))
    rd <- 0
  return(rd)
} 

##------------------------------------------------------------- 
##Separate variables into precip and temperature
## 'E' denotes a climdex index
##Set the order that the variables should be arranged here.

filter.input.variables <- function(table.vars) {
  all.vars <- unlist(lapply(table.vars,function(x){return(x[1])}))
  pr.vars <- c('pr.ann.avg','pr.ann.max',
               'pr_rp5','pr_rp20','pr_rp50',
               'snow_rp20','snow_rp50')
  pr.ix <- match(pr.vars,all.vars)
  pr.selected <- table.vars[pr.ix[!is.na(pr.ix)]]

  tas.vars <- c('tasmax.950.mon','tasmax.975.mon','tasmax.990.mon','tasmax.996.mon',
                'tasmin.004.mon','tasmin.010.mon','tasmin.025.mon','tasmin.050.mon',
                'wb.950.mon','wb.975.mon','wb.990.mon','wb.996.mon',
                'cdd.18.ann','cdd.18.3.ann','hdd.18.ann','hdd.18.3.ann')
  tas.ix <- match(tas.vars,all.vars)
  tas.selected <- table.vars[tas.ix[!is.na(tas.ix)]]

  other.vars <- c('wind_rp5','wind_rp10','wind_rp50')
  other.ix <- match(other.vars,all.vars)
  other.selected <- table.vars[other.ix[!is.na(other.ix)]]
  rv <- list(pr=pr.selected,tas=tas.selected,other=other.selected)

  return(rv)
}

##------------------------------------------------------------- 

find.row.locations <- function(sorted.vars) {

   pr.vars <- sorted.vars$pr
   pr.rows <- vector(length=length(pr.vars)+1,mode='list')

   pr.rows[[1]] <- 3:5 ##Header Rows
   for (i in 1:length(pr.vars)) {
      row.len <- switch(pr.vars[[i]][[2]],annual=2,seasonal=6)
      pr.rows[[i+1]] <- seq(tail(pr.rows[[i]],1)+1,length.out=row.len)
   }

   tas.start <- tail(pr.rows[[length(pr.vars)+1]],1)
   tas.vars <- sorted.vars$tas
   tas.rows <- vector(length=length(tas.vars)+1,mode='list')
   tas.rows[[1]] <- (1:3)+tas.start
   for (i in 1:length(tas.vars)) {
      row.len <- switch(tas.vars[[i]][[2]],annual=2,seasonal=6)
      tas.rows[[i+1]] <- seq(tail(tas.rows[[i]],1)+1,length.out=row.len)
   }

   other.start <- tail(tas.rows[[length(tas.vars)+1]],1)
   other.vars <- sorted.vars$other
   other.rows <- vector(length=length(other.vars)+1,mode='list')
   other.rows[[1]] <- (1:3)+other.start
   for (i in 1:length(other.vars)) {
      row.len <- switch(other.vars[[i]][[2]],annual=2,seasonal=6)
      other.rows[[i+1]] <- seq(tail(other.rows[[i]],1)+1,length.out=row.len)
   }

   rv <- list(pr=pr.rows,tas=tas.rows,other=other.rows)
   return(rv)

}
##---------------------------------------------

annual.table <- function(var.name,scenario='rcp85',rp=NULL) {
  
  no.percent <- '(tas|tasmax|tasmin|txxETCCDI|tnnETCCD|trETCCDI|r95sep|r99days|r95days)'
  result <- get.annual.data(var.name,scenario,rp) 

  if (grepl(no.percent,var.name))
    result[9:14] <- 'NA'
  return(as.data.frame(t(result)))
}


##----------------------------------------------------------------------------------------
write.variables <- function(wb,sorted.vars,row.locs,type) {
  len <- length(sorted.vars)
  for (i in 1:len) {               
    current.var <- sorted.vars[[i]]
    var.name <- current.var[1]
    print(var.name)
    season <- current.var[2]
    var.title <- current.var[3]
    print(var.title)
    current.row <- row.locs[[i+1]]

    s1 <- createStyle(fontSize = 12, fontColour = "black", textDecoration = c("BOLD"))
    s2 <- createStyle(fontSize = 12, fontColour = "black")
    c1 <- createComment(comment = variable.comment(var.name),style=c(s1,s2),visible=FALSE)
    writeComment(wb, 1, col = 1, row = current.row[1], comment = c1)
    seas.fx <- switch(season,   
                      annual=annual.table,   
                      seasonal=seasonal.table)
    if (grepl('rp',var.name)) {
      rp <- as.numeric(gsub('rp','',strsplit(var.name,'_')[[1]][2]))
      rpvar <- gsub('rp','',strsplit(var.name,'_')[[1]][1])
      var.entry <- seas.fx(rpvar,rp=rp)
    } else {                      
      var.entry <- seas.fx(var.name)                      
    }
    pane.colour <- switch(type,pr='lightblue',tas='tan1')             

    hdstyle <- createStyle(fgFill = pane.colour, halign = "CENTER", textDecoration = "Bold",
                             border = "TopBottomLeftRight", fontColour = "black")  
    units.pane <- c(var.title,get.units.pane(var.name))                             
    writeData(wb, sheet=1, as.data.frame(t(units.pane)), startRow = current.row[1], startCol = 1, headerStyle = hdstyle,
              borders = "rows", borderStyle = "medium",colNames=FALSE)
    addStyle(wb,sheet=1,hdstyle,rows=current.row[1],cols=1:14,gridExpand=FALSE,stack=FALSE)

    datastyle <- createStyle(fgFill = 'white', halign = "RIGHT",
                             border = "TopBottomLeftRight", fontColour = "black")                              
    writeData(wb, sheet=1, var.entry, startRow = current.row[2], startCol = 1, headerStyle = hdstyle,
              borders = "rows", borderStyle = "medium",colNames=FALSE)
    addStyle(wb,sheet=1,datastyle,rows=current.row[-1],cols=1:14,gridExpand=TRUE,stack=FALSE)

    highlight <- createStyle(fgFill = 'lightyellow', halign = "CENTER", 
                             border = "TopBottomLeftRight", fontColour = "black")  
    addStyle(wb,sheet=1,highlight,rows=current.row[-1],cols=2,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=1,highlight,rows=current.row[-1],cols=5,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=1,highlight,rows=current.row[-1],cols=6,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=1,highlight,rows=current.row[-1],cols=11,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=1,highlight,rows=current.row[-1],cols=12,gridExpand=FALSE,stack=FALSE)

    if (grepl('(suETCCDI|su30ETCCDI)',var.name) | var.name =='cdd') {
        prctstyle <- createStyle(fgFill = 'lightgray', halign = "RIGHT",
                                 border = "TopBottomLeftRight", fontColour = "black")                              
        writeData(wb, sheet=1, var.entry[9:14], startRow = current.row[2], startCol = 9, headerStyle = prctstyle,
                  borders = "rows", borderStyle = "medium",colNames=FALSE)
        addStyle(wb,sheet=1,prctstyle,rows=current.row[-1],cols=9:14,gridExpand=TRUE,stack=FALSE)
        s1 <- createStyle(fontSize = 12, fontColour = "black", textDecoration = c("BOLD"))
        s2 <- createStyle(fontSize = 12, fontColour = "black")
        f1 <- createComment(comment = c('CAUTION\n','Percent changes from a low baseline value can result in deceptively large percent values'),
                            style=c(s1,s2),visible=FALSE)
        for (j in 9:14) {                     
          writeComment(wb, 1, col = j, row = current.row[1]+1, comment = f1)
        }

    }


  }
     
}
##----------------------------------------------------------------------------------------
      ##Top Frozen Pane

create.frozen.top.pane <- function(wb) {

      pane.titles <- list(' ','Past','2020s Change',' ',
                              '2050s Change',' ',
                              '2080s Change',' ',
                              '2020s Percent Change',' ',
                              '2050s Percent Change',' ',
                              '2080s Percent Change',' ')

      fz1 <- createStyle(fgFill = "gray94", halign = "CENTER", textDecoration = "Bold",
                         border = "TopBottomLeftRight", fontColour = "black")
      writeData(wb, sheet=1, pane.titles, startRow = 1, startCol = 1, headerStyle = fz1,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=1,fz1,rows=1,cols=1:14,gridExpand=FALSE,stack=FALSE)
      mergeCells(wb,sheet=1,cols=c(3,4),rows=1)
      mergeCells(wb,sheet=1,cols=c(5,6),rows=1)
      mergeCells(wb,sheet=1,cols=c(7,8),rows=1)
      mergeCells(wb,sheet=1,cols=c(9,10),rows=1)
      mergeCells(wb,sheet=1,cols=c(11,12),rows=1)
      mergeCells(wb,sheet=1,cols=c(13,14),rows=1)
      ##freezePane(wb,sheet=1,firstRow=TRUE)

      prct.header <- list(' ',' ','Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%',
                                  'Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%')          
      writeData(wb, sheet=1, prct.header, startRow = 2, startCol = 1, headerStyle = fz1,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=1,fz1,rows=2,cols=1:14,gridExpand=FALSE,stack=FALSE)
      ##freezePane(wb,sheet=1,firstActiveRow=3)
}

##-----------------------------------------------------------------------------------------
      ##Precipitation Header Rows
create.title.panes <- function(wb,var.name,start.row) {

      pane.titles <- list(' ','Past','2020s Change',' ',
                               '2050s Change',' ',
                              '2080s Change',' ',
                              '2020s Percent Change',' ',
                              '2050s Percent Change',' ',
                              '2080s Percent Change',' ')

      pane.colour <- switch(var.name,pr='lightblue',tas='tan1')             
      hdstyle <- createStyle(fgFill = pane.colour, halign = "CENTER", textDecoration = "Bold",
                             border = "TopBottomLeftRight", fontColour = "black") 
      writeData(wb, sheet=1, pane.titles, startRow = start.row, startCol = 1, headerStyle = hdstyle,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=1,hdstyle,rows=start.row,cols=1:14,gridExpand=FALSE,stack=FALSE)
      mergeCells(wb,sheet=1,cols=c(3,4),rows=start.row)
      mergeCells(wb,sheet=1,cols=c(5,6),rows=start.row)
      mergeCells(wb,sheet=1,cols=c(7,8),rows=start.row)
      mergeCells(wb,sheet=1,cols=c(9,10),rows=start.row)
      mergeCells(wb,sheet=1,cols=c(11,12),rows=start.row)
      mergeCells(wb,sheet=1,cols=c(13,14),rows=start.row)
      prct.header <- list(' ',' ','Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%',
                                  'Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%')          
      writeData(wb, sheet=1, prct.header, startRow = start.row+1, startCol = 1, headerStyle = hdstyle,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=1,hdstyle,rows=start.row+1,cols=1:14,gridExpand=FALSE,stack=FALSE)

      mergeCells(wb,sheet=1,cols=1:14,rows=start.row+2)
      pane.header <- list(switch(var.name,pr='Precipitation',tas='Temperature'))
      writeData(wb, sheet=1, pane.header, startRow = start.row+2, startCol = 1, headerStyle = hdstyle,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=1,hdstyle,rows=start.row+2,cols=1:14,gridExpand=FALSE,stack=FALSE)
}


##*******************************************************************
##
##*******************************************************************

table.vars <-    list(c('pr_rp5','annual','RP5 PR'),
                      c('pr_rp20','annual','RP20 PR'),
                      c('pr_rp50','annual','RP50 PR'),
                      c('pr.ann.avg','annual','ANN PR AVG'),
                      c('pr.ann.max','annual','ANN PR MAX'),
                      c('snow_rp20','annual','RP20 SNOW'),
                      c('snow_rp50','annual','RP50 SNOW'),
                      c('tasmax_rp5','annual','RP5 TX'),
                      c('tasmin_rp5','annual','RP5 TN'),
                      c('tasmax.950.mon','annual','TX 95.0%'),
                      c('tasmax.975.mon','annual','TX 97.5%'),
                      c('tasmax.990.mon','annual','TX 99.0%'),
                      c('tasmax.996.mon','annual','TX 99.6%'),
                      c('tasmin.004.mon','annual','TN 0.4%'),
                      c('tasmin.010.mon','annual','TN 1.0%'),
                      c('tasmin.025.mon','annual','TN 2.5%'),
                      c('tasmin.050.mon','annual','TN 5.0%'),
                      c('wb.950.mon','annual','WB 95.0%'),
                      c('wb.975.mon','annual','WB 97.5%'),
                      c('wb.990.mon','annual','WB 99.0%'),
                      c('wb.996.mon','annual','WB 99.6%'),
                      c('cdd.18.ann','annual','CDD 18'),
                      c('cdd.18.3.ann','annual','CDD 18.3'),
                      c('hdd.18.ann','annual','HDD 18'),
                      c('hdd.18.3.ann','annual','HDD 18.3'),
                      c('wind_rp5','annual','RP5 WIND'),
                      c('wind_rp10','annual','RP10 WIND'),
                      c('wind_rp50','annual','RP50 WIND'))
                       
##------------------------------------------------------------------------


write.formatted.table <- function(bcbc.entries,bcbc.names,table.vars) {

  sorted.vars <- filter.input.variables(table.vars) 

  row.locs <- find.row.locations(sorted.vars)

  ##Formatted Table
  wb <- createWorkbook()
  addWorksheet(wb, "Regional Averages")
  setColWidths(wb, sheet = 1, cols = 1:14, widths = 14) ##Set fixed width
  create.frozen.top.pane(wb)
  create.title.panes(wb,var.name='pr',start.row=row.locs$pr[[1]][1])
  write.variables(wb,sorted.vars$pr,row.locs$pr,'pr')
  create.title.panes(wb,var.name='tas',start.row=row.locs$tas[[1]][1])
  write.variables(wb,sorted.vars$tas,row.locs$tas,'tas')
  create.title.panes(wb,var.name='other',start.row=row.locs$other[[1]][1])
  write.variables(wb,sorted.vars$other,row.locs$other,'other')


  freezePane(wb,sheet=1,firstActiveCol=2,firstActiveRow=3)

  ##saveWorkbook(wb, paste0('/storage/home/ssobie/general/assessment_tables/',region,'_variable_table_rcp85.xlsx'), overwrite = TRUE)
  saveWorkbook(wb, paste0('/storage/data/projects/rci/building_code/',readloc,'/tables/',region,'_variable_table_rcp85.xlsx'), overwrite = TRUE)

}

