####################
## PMF execution and postprocessing program
## ~userinputs.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################

## CalNex

FOLDER <- "runs/FTIR02"
nFactors <- 12##c(2,5,8)
FPEAK <- round(seq(-1.2,.6,1.2),2)
Seeds <- c(1,100)
newsim <- FALSE

