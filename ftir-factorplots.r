####################
## PMF execution and postprocessing program
## ~postprocess.r~
## $Rev: -1 $
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


## input:
input <- commandArgs()
pattern <- "--file=(.+)"
srcpath <- gsub('~+~'," ",dirname(sub(pattern,"\\1",input[grepl(pattern,input)])),fixed=TRUE)
source(file.path(srcpath,"functions/io.R"))

argv <- tail(input,-grep("--args",input,fixed=TRUE))
filename <- argv[1]

args <- readArgs(filename)
for(p in names(args))
  assign(p,args[[p]])

source(file.path(srcpath,"functions/ftir-factorplots.r"))
patt <- ".+\\_([0-9]{3})\\.pdf"
######################################################
# output: Allplots/Allplots.pdf
######################################################
setwd(FOLDER)

# post-process
dir.create("Allplots")
## source(file.path(homedir,"functions/results_template.r"))
## lapply(c("Qvalues.pdf","Gmedian.pdf"),
##        function(x) file.rename(x,file.path("Allplots",x)))

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

