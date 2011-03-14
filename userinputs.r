####################
## PMF execution and postprocessing program
## ~userinputs.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


## CalNex
FOLDER <- "runs/ACSM02"
nFactors <- c(2,8)
FPEAK <- round(seq(-1.2,1.2,1.2),2)
Seeds <- c(1,100)
newsim <- TRUE
xvariables <- "amus.txt"

## FOLDER <- "runs/test"
## nFactors <- 2
## FPEAK <- 1
## Seeds <- 1
## newsim <- TRUE
