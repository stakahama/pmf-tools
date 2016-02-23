

FOLDER <- "~/VirtualBoxDocuments/programs/pmf"

addstdev <- function(run,runno,newrun,add) {
  ##
  runno <- sprintf("%s_%03d",basename(run),runno)
  simgrid <- read.delim(file.path(FOLDER,run,"simgrid.txt"),row.names=1)
  ##
  ff <- as.matrix(read.table(file.path(FOLDER,run,runno,"F_FACTOR.TXT")))
  gg <- as.matrix(read.table(file.path(FOLDER,run,runno,"G_FACTOR.TXT")))
  noise <- do.call(`+`,lapply(add,function(i,gg,ff) gg[,i] %*% t(ff[i,]),gg,ff))
  ##
  stdev <- as.matrix(read.table(file.path(FOLDER,run,"STD_DEV.DAT")))
  ##
  dir.create(file.path(FOLDER,newrun))
  write.table(stdev+noise,file.path(FOLDER,newrun,"STD_DEV.DAT"),
              row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
  for( x in c("MATRIX.DAT","samples.txt","variables.txt") )
    file.copy(file.path(FOLDER,run,x),file.path(FOLDER,newrun,x))

}

addstdev(run="runs/ACSM10",
         runno=16,
         newrun="runs/ACSM16",
         2:3)

addstdev(run="runs/ACSM16",
         runno=2,
         newrun="runs/ACSM17",
         add=3)
