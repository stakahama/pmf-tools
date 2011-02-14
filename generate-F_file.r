####################
## PMF execution and postprocessing program
## ~generate-F_file.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


source("userinputs.r")
ffiles <- list.files(FOLDER,"F_FACTOR.TXT",recursive=TRUE,full=TRUE)

wavenums <- scan(file.path(FOLDER,"wavenumbers.txt"),0,quiet=TRUE)
outfile <- "runs/FIN3to6/allfactors.txt"

cat(paste(c("",formatC(wavenums,format="f",digits=3)),collapse="\t"),"\n",
    sep="",file=outfile)
for( x in ffiles ) {
  fmat <- matrix(scan(x,0,quiet=TRUE),ncol=length(wavenums),byrow=TRUE)
  if( nrow(fmat) == 0 ) next
  rownames(fmat) <- sprintf("%s-%02d",basename(dirname(x)),1:nrow(fmat))
  write.table(fmat,file=outfile,sep="\t",col.names=FALSE,
              quote=FALSE,append=TRUE)
}

