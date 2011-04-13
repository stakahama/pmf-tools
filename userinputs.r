####################
## PMF execution and postprocessing program
## ~userinputs.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################

## CalMex

FOLDER <- "runs/ACSM07"
nFactors <- 3:4
FPEAK <- round(seq(-1.2,1.2,1.2),2)
Seeds <- c(1,100)
newsim <- FALSE

## FOLDER <- "runs/FTIR02"
## nFactors <- c(2,3,6,8)
## FPEAK <- 1.2
## Seeds <- c(1,10,100)
## newsim <- FALSE
