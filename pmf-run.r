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
simgrid <- local({
  new <- expand.grid(nFactors=nFactors,FPEAK=FPEAK,Seed=Seeds)
  new <- new[with(new,order(nFactors,FPEAK)),]
  patt <- sprintf("%s\\_([0-9]+)",basename(FOLDER))
  mx <- (if( newsim ) 0 else
         (if(file.exists(FOLDER) &&
             length(fi <- list.files(FOLDER,patt,full=TRUE)) > 0)
          max(0,as.integer(na.omit(sub(patt,"\\1",
                                       basename(fi[file.info(fi)$isdir])))),
              na.rm=TRUE) else 0))
  `rownames<-`(new,paste(basename(FOLDER),
                         substring(1000+mx+(1:nrow(new)),2),sep="_"))
})

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
