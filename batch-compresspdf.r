
source("userinputs.r")

compresspdf <- function(x) {
  ## requires that pdftk is installed
  tmp <- paste(x,"~",sep="")
  system(sprintf("pdftk %s output %s compress dont_ask",x,tmp))
  file.rename(tmp,x)
  cat(paste("compressed",x),"\n")  
}

pdflist <- list.files(FOLDER,"\\.pdf",rec=TRUE,full=TRUE)
invisible(Map(compresspdf,pdflist))
