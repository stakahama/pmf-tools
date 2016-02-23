postprocess <- function() {
  source("pmf-qvalues.r",local=TRUE)
  source("pmf-mutual-correlations.r",local=TRUE)    
  source("pmf-explvar.r",local=TRUE)      
  source("pmf-explvar-plots.r",local=TRUE)
}

## inputs <- readLines("userinputs.r")
## runs <- list.files("runs","532\\_",full=TRUE)

## for( x in runs ) {
##   print(x)
##   newlines <- replace(inputs,grep("^(?!## )FOLDER",inputs,perl=TRUE),
##                       sprintf('FOLDER <- "%s"',x))
##   writeLines(newlines,"userinputs.r")
##   tryCatch({
##     postprocess()
##     ## source("stxm-factorplots.r",local=TRUE)    
##   },error=function(e) print(e))
## }

postprocess()
## system("Rscript pmf-g-correlations.r")
## system("Rscript acsm-factor-classification.r")
system("Rscript pmf-g-correlations.r")
system("Rscript ftir-factor-classification.r")        
