####################
## PMF execution and postprocessing program
## ~g-correlations.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


###_* user inputs
input <- commandArgs()
pattern <- "--file=(.+)"
srcpath <- gsub('~+~'," ",dirname(sub(pattern,"\\1",input[grepl(pattern,input)])),fixed=TRUE)
source(file.path(srcpath,"functions/io.R"))

argv <- tail(input,-grep("--args",input,fixed=TRUE))
filename <- argv[1]

###_* libraries
library(lattice)
library(latticeExtra)

args <- read.args(filename)
for(p in names(args))
  assign(p,args[[p]])

###_* functions
read.G <- function(FOLDER,run) {
  g <- as.matrix(read.table(file.path(FOLDER,run,"G_FACTOR.TXT")))
  dimnames(g) <- list(readLines(file.path(FOLDER,"samples.txt")),
                      sprintf("%s-%02d",run,1:ncol(g)))
  g
}
icor <- function(run) {
  x <- cor(read.G(FOLDER,run))
  y <- outer(rownames(x),colnames(x),paste,sep=":")
  data.frame(run=run,pair=y[upper.tri(y)],r=x[upper.tri(x)])
}

###_* simgrid
simgrid <- read.delim(file.path(FOLDER,"simgrid.txt"),row.names=1)
simgrid <- simgrid[with(simgrid,order(nFactors,FPEAK,Seed)),]
## sims <- rownames(subset(simgrid,nFactors %in% 2:6))
sims <- rownames(simgrid)

###_* calculations
correl <- do.call(rbind,lapply(sims,icor))
dfr <- cbind(correl,simgrid[correl$run,])
dfr <- within(dfr,{
  Seed <- factor(Seed)
  nFactors <- factor(nFactors)
  ## FPEAK <- factor(FPEAK)
})

###_* plot and export

panel <- function(...) {
  panel.abline(h=seq(-1,1,.5),col=8,lty=2)  
  panel.xyplot(...)    
}
out <- xyplot(r~FPEAK|nFactors*Seed,data=dfr,
              scales=list(x=list(rot=90),y=list(at=seq(-1,1,.5))),
              panel=panel,cex=.5,
              as.table=TRUE,ylim=c(-1,1))
pdf(file.path(FOLDER,"Allplots","G-correlations.pdf"),width=8,height=5)
useOuterStrips(out)
dev.off()
