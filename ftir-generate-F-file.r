####################
## PMF execution and postprocessing program
## ~generate-F_file.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


source("userinputs.r")
## get all files
ffiles <- list.files(FOLDER,"F_FACTOR.TXT",recursive=TRUE,full=TRUE)

## read in wavenumbers
wavenums <- scan(file.path(FOLDER,"variables.txt"),0,quiet=TRUE)

## write allfiles
for( x in ffiles ) {
  fmat <- matrix(scan(x,0,quiet=TRUE),ncol=length(wavenums),byrow=TRUE)
  if( nrow(fmat) == 0 ) next
  outfile <- file.path(dirname(x),"allfactors.txt")  
  rownames(fmat) <- sprintf("%s-%02d",basename(dirname(x)),1:nrow(fmat))
  cat(paste(c("",formatC(wavenums,format="f",digits=3)),collapse="\t"),"\n",
      sep="",file=outfile)
  write.table(fmat,file=outfile,sep="\t",col.names=FALSE,
              quote=FALSE,append=TRUE)
}

