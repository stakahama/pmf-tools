####################
## PMF execution and postprocessing program
## ~userinputs.r~
## $Rev: 9 $
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################

## CalNex

FOLDER <- "runs/ACSM07"
nFactors <- 20
FPEAK <- round(seq(-1.2,1.2,1.2),2)
Seeds <- c(1,100)
newsim <- FALSE
