####################
## PMF execution and postprocessing program
## ~exportt2dpt.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


## User inputs
source("userinputs.r")
dptpath <- "dpt"

## Define DPT generating function
write.dpt <- function(path,FOLDER) {
  ## wavenums and dptpath are global
  F <- matrix(scan(file.path(FOLDER,path,"F_FACTOR.TXT"),0,quiet=TRUE),
              ncol=length(wavenums),byrow=TRUE)
  ifactor <- 1
  for( f in split(F,row(F)) ) {
    write(rbind(wavenums,f), ## intentionally transposed
          file=file.path(dptpath,sprintf("%s-%02d.0.DPT",basename(path),ifactor)),
          ncolumns=2,sep=",")
    ifactor <- ifactor + 1
  }
}

## scan wavenumbers, create DPT path
wavenums <- scan(file.path(FOLDER,"wavenumbers.txt"),0,quiet=TRUE)
dir.create(dptpath)

## Apply function
runs <- list.files(FOLDER,basename(FOLDER))
invisible(lapply(runs[1:2],write.dpt,FOLDER))
