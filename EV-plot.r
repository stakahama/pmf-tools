####################
## PMF execution and postprocessing program
## ~EV-plot.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


ev <- read.delim("outputs/ExplainedVariance-TJ1_012.txt",check.names=FALSE)
ev <- read.delim("outputs/ExplainedVariance-TJ1_016.txt",check.names=FALSE)

norm <- function(x) x/sum(x)
ev[] <- t(apply(ev,1,norm))

library(fields)
barplot(t(ev),legend=TRUE,col=tim.colors(ncol(ev)))
