####################
## PMF execution and postprocessing program
## ~run_template.r~
## $Rev: -1 $
## Feb. 2016
## Satoshi Takahama (satoshi.takahama@epfl.ch)
####################

## inputs:
input <- commandArgs()
pattern <- "--file=(.+)"
srcpath <- gsub('~+~'," ",dirname(sub(pattern,"\\1",input[grepl(pattern,input)])),fixed=TRUE)
source(file.path(srcpath,"functions/io.R"))

argv <- tail(input,-grep("--args",input,fixed=TRUE))
filename <- argv[1]

## contents
args <- readArgs(filename)

#######################################################
# output : lots of files
#######################################################
source(file.path(srcpath,"functions/Rfuncs.r"))

## construct $X$, $\sigma$ ** (many options)
X <- as.matrix(read.table(file.path(args[["FOLDER"]],"matrix.dat")))
Xstdev <- as.matrix(read.table(file.path(args[["FOLDER"]],"std_dev.dat")))

## set up grid
if(is.null(args[["newsim"]]))
  args[["newsim"]] <- TRUE
simgrid <- with(args,creategrid(nFactors,FPEAK,Seed,FOLDER,newsim))

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
    if(!file.exists(args[["FOLDER"]])) dir.create(args[["FOLDER"]])
    try(setwd(args[["FOLDER"]]),TRUE)
  }
  write.table(simgrid,file="simgrid.txt",sep="\t",col=NA,quote=FALSE,
              append=if(args[["newsim"]]) FALSE else TRUE)  

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
