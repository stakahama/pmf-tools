####################
## PMF execution and postprocessing program
## ~postprocess.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


## input:
source("userinputs.r")
patt <- ".+\\_([0-9]{3})\\.pdf"
######################################################
# output: Allplots/Allplots.pdf
######################################################
homedir <- getwd()
setwd(FOLDER)

# post-process
try(dir.create("Allplots"),TRUE)
## source(file.path(homedir,"functions/results_template.r"))
## lapply(c("Qvalues.pdf","Gmedian.pdf"),
##        function(x) file.rename(x,file.path("Allplots",x)))
source(file.path(homedir,"functions/ftir-factorplots.r"))
setwd("Allplots")
try({
  projpdf <- list.files(".",patt=patt)
  pdfSeq <- c(##"Qvalues.pdf","Gmedian.pdf",
              projpdf[order(as.numeric(gsub(patt,"\\1",projpdf)))])
  outName <- "Factorplots.pdf"
  file.remove(outName)
  cmd <- c("pdftk",pdfSeq,"cat output",outName,"compress")
  output <- system(paste(cmd,collapse=" "),intern=TRUE)
  lapply(pdfSeq,file.remove)
})
setwd(homedir)
