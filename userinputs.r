####################
## PMF execution and postprocessing program
## ~userinputs.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


## CalNex
FOLDER <- "runs/ACSM01"
nFactors <- 2:8
FPEAK <- seq(-1.2,1.2,.4)
Seeds <- c(1,10,100)
newsim <- TRUE
xvariables <- "amus.txt"

## FOLDER <- "runs/test"
## nFactors <- 2
## FPEAK <- 1
## Seeds <- 1
## newsim <- TRUE
