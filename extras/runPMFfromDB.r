####################
## PMF execution and postprocessing program
## ~runPMFfromDB.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


library(RSQLite)
dbname <- "../clusterspec/particledb.db"  
drv <- dbDriver("SQLite")
con <- dbConnect(drv,dbname=dbname)
tblnms <- dbGetQuery(con,"SELECT * from SQLITE_MASTER")$tbl_name
tblnms <- c('P60716205')
cdir <- getwd()
for( tb in tblnms ) {
  FOLDER <- sprintf("%sp",substring(tb,2))
  if( file.exists(FOLDER) ) next
  tryCatch({
    X <-
      dbGetQuery(con,sprintf('SELECT %s FROM %s WHERE bw = 1',
                             paste(sprintf('en%d',1:95),collapse=','),tb))
##     X[] <- lapply(X,function(x)
##                   if(is.character(x))
##                   type.convert(x,na.string='"NULL"',as.is=TRUE)
##                   else x)
    Xstdev <- replace(X,TRUE,1)
    dir.create(FOLDER)
    write.table(round(X,8),file=file.path(FOLDER,"matrix.dat"),
                sep="\t",col=FALSE,row=FALSE)
    write.table(round(Xstdev,8),file=file.path(FOLDER,"std_dev.dat"),
                sep="\t",col=FALSE,row=FALSE)
    source("run_template.r")
    source("postprocess.r")
  },error=function(e) NULL)
  setwd(cdir)
}
dbDisconnect(con)


