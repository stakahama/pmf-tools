####################
## PMF execution and postprocessing program
## ~userinputs.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################

## CalMex

## FOLDER <- "runs/ACSM12"
## nFactors <- 2
## FPEAK <- c(-.8,-.4)
## Seeds <- c(1,10,100)

## FOLDER <- "runs/ACSM09"
## nFactors <- 3
## FPEAK <- 0
## Seeds <- 1

## FOLDER <- "runs/FTIR03"
## nFactors <- c(3:5)
## FPEAK <- round(seq(-1.2,1.2,.6),2)
## Seeds <- c(1,10,100)
## newsim <- TRUE

FOLDER <- "runs/FTIR04"
nFactors <- 4
FPEAK <- c(-1.8,1.8)
Seeds <- c(1,10,100)

## FOLDER <- "runs/532_091125014_p1"
## nFactors <- 2:3
## FPEAK <- c(-1.2,-0.6,0,0.6,1.2)
## Seeds <- c(1,100)
