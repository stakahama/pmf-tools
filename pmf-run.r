####################
## PMF execution and postprocessing program
## ~run_template.r~
## $Rev: -1 $
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


## inputs:
source("userinputs.r")
#######################################################
# output : lots of files
#######################################################
source("functions/Rfuncs.r")

## construct $X$, $\sigma$ ** (many options)
X <- as.matrix(read.table(file.path(FOLDER,"matrix.dat")))
Xstdev <- as.matrix(read.table(file.path(FOLDER,"std_dev.dat")))

## set up grid
newsim <- (if(file.exists(file.path(FOLDER,"simgrid.txt"))) FALSE else TRUE)
simgrid <- creategrid(nFactors,FPEAK,Seeds,newsim)

## (this is a function which requires nFactors and FPEAK)
## also make sure $X$ is created and exists in global space
jointrepl <-
    with(list(lines=readLines("c:/PMF/pmf2def.ini")),makeFunc(ini=lines,mat=X))
##with(list(lines=readLines("pmf2def.ini")),makeFunc(ini=lines,mat=X))

## write $X$, $\sigma$, and simgrid
homedir <- getwd()
if( basename(getwd()) != FOLDER ) {
  try(dir.create(FOLDER),TRUE)
  try(setwd(FOLDER),TRUE)
}
if(newsim) {
  write.table(simgrid,file="simgrid.txt",sep="\t",col=NA,quote=FALSE)  
} else {
  write.table(simgrid,file="simgrid.txt",sep="\t",col=FALSE,quote=FALSE,
              append=TRUE)
}

##run 
invisible(sapply(split(simgrid,simgrid$nFactors),
                 function(x) sapply(1:nrow(x),execPMF,x)))
lapply(c("F_FACTOR.TXT","G_FACTOR.TXT","TEMP.TXT","MISC.TXT","PMF2.LOG",
         "MYPMF.INI"),file.remove)

setwd(homedir)
