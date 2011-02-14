####################
## PMF execution and postprocessing program
## ~qvalues.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


###_* read simgrid

###_ . whistler
source("userinputs.r")
simgrid <- read.delim(file.path(FOLDER,"simgrid.txt"),row.names=1,as.is=TRUE)

###_* read dim(X), Q-values

###_ . X
Xmat <- read.table(file.path(FOLDER,"matrix.dat"))
eqm <- function(p) prod(dim(Xmat)) - p*(sum(dim(Xmat)))

###_ . Q-values
outputs <- list.files(FOLDER,"output.txt",recursive=TRUE,full=TRUE)
maybe <- function(x) if(length(x)==0) NA else x
qval <- function(x, patt=".*Final chi2[=][ ]*([0-9.]+).*")
  maybe(as.numeric(sub(patt,"\\1",grep(patt,readLines(x),value=TRUE))))
vals <- sapply(outputs,qval)
names(vals) <- basename(dirname(names(vals)))

###_* plot
library(lattice)
simgrid$`Q-ratio` <- vals[rownames(simgrid)]/eqm(simgrid$nFactors)
simgrid$nFactors <- factor(simgrid$nFactors)

panel <- function(...) {
  panel.abline(h=0,col=8)
  panel.xyplot(...)
}
xout <- xyplot(`Q-ratio`~FPEAK | nFactors, groups=Seed,data=simgrid,
               panel=panel,scales=list(y=list(relation="sliced")),
               type="o",auto=list(space="right",title="Seed",cex.title=1),
               ylab=expression(Q/Q[expected]),
               as.table=TRUE)

pdf(file.path(FOLDER,"Allplots","Qvalues.pdf"),width=8,height=5)
## print(update(xout,scales=list(y=list(relation="same"))))
print(xout)
dev.off()
