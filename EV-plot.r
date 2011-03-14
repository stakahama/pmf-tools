####################
## PMF execution and postprocessing program
## ~EV-plot.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################

###_* user inputs
source("userinputs.r")

###_* single

## norm <- function(x) x/sum(x)
## ev <- read.delim(file.path(FOLDER,"ACSM01_076","ExplainedVariation.txt"),check.names=FALSE)
## ev[] <- t(apply(ev,1,norm))
## library(fields)
## barplot(t(ev),legend=TRUE,col=tim.colors(ncol(ev)),border=NA)

###_* multiple

###_ . explained variation statistics

sumstats <- function(filename) {
  ev <- read.delim(filename,row.names=1)
  x <- 1-ev$Resid
  c(mean=mean(x),median=median(x),min=min(x),max=max(x))
}
evstats <- t(sapply(list.files(FOLDER,"ExplainedVariation.txt",rec=TRUE,full=TRUE),sumstats))
rownames(evstats) <- basename(dirname(rownames(evstats)))

###_ . merge with simgrid
simgrid <- read.delim(file.path(FOLDER,"simgrid.txt"),row.names=1)
simgrid <- with(list(x=merge(simgrid,data.frame(evstats),by="row.names")),
                `row.names<-`(x[,-1],x[,1]))
simgrid$Seed <- factor(simgrid$Seed)
simgrid$FPEAK <- factor(round(simgrid$FPEAK,2))

###_* plot

library(lattice)
tpars <- trellis.par.get("superpose.symbol")

pdf(file.path(FOLDER,"Allplots","ExplainedVariation-nFactors.pdf"),width=8,height=5)
par(mar=c(4,4,1,5),mgp=c(2.2,.5,.0))
with(simgrid,{
  plot(nFactors,median*100,col=tpars$col[Seed],
       pch=unclass(FPEAK),
       ylim=c(0,100),ylab="Explained Variation (%)")
  legend(par("usr")[2],par("usr")[4],xjust=0,yjust=1,
         title="Seed",pch=1,col=tpars$col[1:nlevels(Seed)],
         levels(Seed),xpd=NA)
  legend(par("usr")[2],mean(par("usr")[3:4]),xjust=0,yjust=1,
         title="FPEAK",pch=1:nlevels(FPEAK),col=1,
         levels(FPEAK),xpd=NA)
})
## with(simgrid,arrows(nFactors,min*100,nFactors,max*100,
##                     angle=90,col=8,length=.01,code=3))
dev.off()
