####################
## PMF execution and postprocessing program
## ~run_template.r~
## $Rev: -1 $
## Feb. 2016
## Satoshi Takahama (satoshi.takahama@epfl.ch)
####################

library(RJSIONIO)

## inputs:
filename <- commandArgs(TRUE)
args <- fromJSON(filename)
if(is.null(args[["FOLDER"]]))
  args[["FOLDER"]] <- path.expand(dirname(filename))

#######################################################
# output : lots of files
#######################################################
source("functions/Rfuncs.r")

## construct $X$, $\sigma$ ** (many options)
X <- as.matrix(read.table(file.path(args[["FOLDER"]],"matrix.dat")))
Xstdev <- as.matrix(read.table(file.path(args[["FOLDER"]],"std_dev.dat")))

## set up grid
newsim <- (if(file.exists(file.path(FOLDER,"simgrid.txt"))) FALSE else TRUE)
simgrid <- with(args,creategrid(nFactors,FPEAK,Seeds,newsim))

## (this is a function which requires nFactors and FPEAK)
## also make sure $X$ is created and exists in global space
jointrepl <-
    with(list(lines=readLines("c:/PMF/pmf2def.ini")),makeFunc(ini=lines,mat=X))
##with(list(lines=readLines("pmf2def.ini")),makeFunc(ini=lines,mat=X))

## main part
homedir <- getwd()
{

  ## write $X$, $\sigma$, and simgrid
  if( basename(getwd()) != args[["FOLDER"]] ) {
    try(dir.create(args[["FOLDER"]]),TRUE)
    try(setwd(args[["FOLDER"]]),TRUE)
  }
  if(newsim) {
    write.table(simgrid,file="simgrid.txt",sep="\t",col=NA,quote=FALSE)  
  } else {
    write.table(simgrid,file="simgrid.txt",sep="\t",col=FALSE,quote=FALSE,
                append=TRUE)
  }

  ##run
  for(x in split(simgrid,simgrid$nFactors)) {
    for(i in 1:nrow(x)) {
      execPMF(i,x)
    }
  }

  ## clean
  lapply(c("F_FACTOR.TXT","G_FACTOR.TXT","TEMP.TXT","MISC.TXT","PMF2.LOG",
           "MYPMF.INI"),file.remove)

}
setwd(homedir)
