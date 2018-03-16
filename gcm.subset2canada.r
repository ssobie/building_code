###Script to assemble a human readable inventory of CMIP5 data at PCIC

slice.by.time <- function(space.file,time.file,gcm,era) {
  if (era) {
    time.write <- gsub(pattern='[0-9]{8}-[0-9]{8}',replacement='19500101-20051231',time.file)
    work <- paste('cdo -O seldate,1950-01-01T00:00,2005-12-31T23:59 ',space.file,' ',time.write,sep='')
  } else {
    time.write <- gsub(pattern='[0-9]{8}-[0-9]{8}',replacement='20060101-21001231',time.file)
    work <- paste('cdo -O seldate,2006-01-01T00:00,2100-12-31T23:59 ',space.file,' ',time.write,sep='')
  }
  system(work)
}

subset.to.canada <- function(input.file,output.file) {
  work <- paste('ncks -O -d lon,215.,310. -d lat,35.,89. ',input.file,' ',output.file,sep='')
  system(work) 
}

concat.files <- function(write.dir,gcm) {
  time.dir <- paste(write.dir,'timetmp/',sep='')
  
  file.list <- list.files(path=time.dir,pattern='*nc')
  hist.file <- file.list[grep('historical',file.list)]
  rcp.files <- file.list[grep('rcp',file.list)]
  
  for (i in seq_along(rcp.files)) {
    gcm.dir <- paste(write.dir,gcm,'/',sep='')
    if (!file.exists(gcm.dir))
      dir.create(gcm.dir,recursive=TRUE)
    rcp.file <- rcp.files[i]
    rst <- regexpr(pattern='rcp[0-9]{2}',rcp.file)
    ren <- rst + attr(rst, "match.length")-1      
    rcp <- substr(rcp.file,rst,ren)
    
    cat.file <- gsub(pattern='[0-9]{8}-[0-9]{8}',replacement='19500101-21001231',hist.file)      
    cat.file <- gsub(pattern='historical',replacement=paste('historical+',rcp,sep=''),cat.file)
    cat.file <- gsub(pattern='time_subset_',replacement='',cat.file)
    work <- paste('ionice -c 3 ncrcat -O ',time.dir,hist.file,' ',time.dir,rcp.file,' ',gcm.dir,cat.file,sep='')
    system(work)

  }
}
  

make.subsets <- function(write.dir) {

  files <- list.files(path=paste(write.dir,'grouptmp',sep=''))
  len <- length(files)
browser()
  for (i in 1:len) {
    file <- files[i]
    era <- grepl('historical',file)
    input.file <- paste(write.dir,'grouptmp/',file,sep='')
    space.file <- paste(write.dir,'spacetmp/space_subset_',file,sep='')
    time.file <- paste(write.dir,'timetmp/time_subset_',file,sep='')
    subset.to.canada(input.file,space.file)
    slice.by.time(space.file,time.file,gcm,era)   
    ##slice.by.time(input.file,time.file,gcm,era)   
  }
}

make.new.files <- function(nc.files,h.ix) {

      ist <- regexpr('[0-9]{8}-',head(nc.files[h.ix],1))
      zst <- ist + attr(ist, "match.length")-2      
      ien <- regexpr('-[0-9]{8}',tail(nc.files[h.ix],1))      
      zen <- ien  + attr(ien, "match.length")-1 

      yst <- substr(head(nc.files[h.ix],1),ist,zst)
      yen <- substr(tail(nc.files[h.ix],1),ien+1,zen)
      new.file <- gsub('[0-9]{8}-[0-9]{8}',paste(yst,'-',yen,sep=''),head(nc.files[h.ix],1))
      return(new.file)
}

group.files <- function(file.list,matrix.subset,write.dir,base.dir,scenarios) {

   ##All scenarios including historical
   for (scenario in scenarios) {
     h.ix <- grep(scenario,file.list)
     if (length(h.ix) > 1) {
        file.paste <- paste(paste(base.dir,file.list[h.ix],sep=''),collapse=' ')
        new.file <- make.new.files(file.list,h.ix)
        new.full <- paste(write.dir,'grouptmp/',new.file,sep='')
        work <- paste('ncrcat -O ',file.paste,' ',new.full,sep='')      

        print(work)
        system(work)
     } else if (length(h.ix)==1) {
        print('File exists')
        work <- paste('cp ',base.dir,file.list[h.ix],' ',write.dir,'grouptmp/',file.list[h.ix],sep='')
        print(work)
        system(work)
     } else if (length(h.ix)==0) {
        print(paste0('No files to concatenate for ',scenario))

     }
   }
}


base.dir <- '/storage/data/climate/downscale/CMIP5/incoming/'

##Make sure the GCM and Centre names are paired in the same order below
gcms <- 'MPI-ESM-MR'

variable.list <- c('pr','tasmax','tasmin')
run <- 'r1i1p1'

full.list <- c('Model','Scenario','Run','Variable','Start','End')
scen.list <- c('historical','rcp26','rcp45','rcp85')

##write.dir <- '/storage/data/climate/downscale/CMIP5/building_code/'

##write.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/CMIP5/global/'
write.dir <- '/storage/data/climate/downscale/BCCAQ2/CMIP5/'

for (gcm in gcms) {
  for (variable in variable.list) {
    ##PR files
    pr.files <- list.files(path=paste(base.dir,gcm,'/download',sep=''),pattern=paste(variable,'_day',sep=''),recursive=TRUE)
    plen <- length(pr.files)
    if (1==1) { ##(plen!=0) {
       ##print(pr.files)
       years <- unlist(regmatches(pr.files,gregexpr('[0-9]{8}-[0-9]{8}',pr.files)))
       yst <- substr(years,1,4)
       yen <- substr(years,10,13)
       
       runs <- unlist(regmatches(pr.files,gregexpr('r[0-9]i1p1',pr.files)))
       rcps <- unlist(regmatches(pr.files,gregexpr('rcp[0-9]{2}|historical',pr.files)))

       file.matrix <- cbind(rep(gcm,plen),rep(variable,plen),rcps,runs,yst,yen)
       test.match <- (file.matrix[,3] %in% scen.list) & (file.matrix[,4] %in% run) 
           
       file.select <- pr.files[test.match]
       matrix.subset <- file.matrix[test.match,]
       print(file.select)
       gcm.dir <- paste0(base.dir,gcm,'/download/')
       ##group.files(file.select,matrix.subset,write.dir,gcm.dir,scen.list)    
   
       make.subsets(write.dir)

       concat.files(write.dir,gcm)
       print('Done concat')

       clean.group <- paste('rm ', write.dir,'grouptmp/*nc',sep='')
       system(clean.group)
       clean.space <- paste('rm ', write.dir,'spacetmp/space_subset*nc',sep='')
       system(clean.space)
       clean.time  <- paste('rm ', write.dir,'timetmp/time_subset*nc',sep='')
       system(clean.time)

    }##If statement
  }##Variable Loop
}##GCM Loop

