
inputs <- readLines("userinputs.r")
runs <- list.files("runs","532\\_",full=TRUE)

for( x in runs ) {
  print(x)
  newlines <- replace(inputs,grep("^(?!## )FOLDER",inputs,perl=TRUE),
                      sprintf('FOLDER <- "%s"',x))
  writeLines(newlines,"userinputs.r")
  tryCatch({
    ## source("qvalues.r",local=TRUE)
    ## source("g-correlations.r",local=TRUE)
    ## source("stxm-cluster-analysis.r",local=TRUE)    
    ## source("EV.r",local=TRUE)      
    ## source("EV-plot.r",local=TRUE)
    source("stxm-factorplots.r",local=TRUE)    
  },error=function(e) print(e))
}
