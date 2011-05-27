####################
## PMF execution and postprocessing program
## ~userinputs.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################

## CalMex

## FOLDER <- "runs/ACSM10"
## nFactors <- 4#2:4
## FPEAK <- -0.6#round(seq(-1.2,1.2,1.2),2)
## Seeds <- c(1,100)
## newsim <- FALSE
## dbname <- "/Users/stakahama/Documents/Work/UCSD/projects/CalMex/programs/acsm-class/allspec_2011-03-12.db"

## FOLDER <- "runs/ACSM09"
## nFactors <- 3
## FPEAK <- 0
## Seeds <- 1
## newsim <- TRUE

FOLDER <- "runs/FTIR03"
nFactors <- c(3:5)
FPEAK <- round(seq(-1.2,1.2,.6),2)
Seeds <- c(1,10,100)
newsim <- TRUE
dbpath <- "/Users/stakahama/Documents/Work/UCSD/projects/CalMex/programs/ftir-class/dbfiles/ftir-refspec.db"
runno <- 18

## FOLDER <- "runs/532_091125014_p1"
## nFactors <- 2:3
## FPEAK <- c(-1.2,-0.6,0,0.6,1.2)
## Seeds <- c(1,100)
## newsim <- TRUE
