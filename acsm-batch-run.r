postprocess <- function() {
  source("pmf-qvalues.r",local=TRUE)
  source("pmf-mutual-correlations.r",local=TRUE)    
  source("pmf-explvar.r",local=TRUE)      
  source("pmf-explvar-plots.r",local=TRUE)
}

postprocess()
system("Rscript pmf-g-correlations.r")
system("Rscript acsm-factor-classification.r")
