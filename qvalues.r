####################
## PMF execution and postprocessing program
## ~qvalues.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


###_* read simgrid

###_ . inputs
source("userinputs.r")
simgrid <- read.delim(file.path(FOLDER,"simgrid.txt"),row.names=1,as.is=TRUE)
dir.create(file.path(FOLDER,"Allplots"))

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
for( var in c("nFactors","Seed") )
  simgrid[,var] <- factor(simgrid[,var])
simgrid$FPEAK <- factor(round(simgrid$FPEAK,2))

panel <- function(...) {
  panel.abline(h=0,col=8)
  panel.xyplot(...)
}
xout <- xyplot(`Q-ratio`~FPEAK | nFactors, groups=Seed,data=simgrid,
               panel=panel,scales=list(y=list(relation="sliced")),
               type="o",auto=list(space="right",title="Seed",cex.title=1),
               ylab=expression(Q/Q[expected]),
               as.table=TRUE)

pdf(file.path(FOLDER,"Allplots","Qvalues-FPEAK.pdf"),width=8,height=5)
## print(update(xout,scales=list(y=list(relation="same"))))
print(xout)
dev.off()

panel <- function(x,y,...,pch,groups,subscripts) {
  panel.points(x,y,pch=pch[subscripts],
               col=trellis.par.get("superpose.symbol")$col[groups[subscripts]])
  panel.lines(smooth.spline(x,y))
}

mykey <- with(simgrid,{
  p <- trellis.par.get("superpose.symbol")
  dfr <- rbind(data.frame(lab=paste("Seed =",levels(Seed)),col=p$col[1:nlevels(Seed)],pch=1),
               data.frame(lab=paste("FPEAK =",levels(FPEAK)),col=1,pch=1:nlevels(FPEAK)))
  list(space="right",text=dfr["lab"],points=dfr[c("col","pch")])
})

xout <- xyplot(`Q-ratio`~nFactors, groups=Seed,data=simgrid,
               panel=panel,pch=simgrid$FPEAK,
               scales=list(y=list(relation="sliced")),
               key=mykey,
               ylab=expression(Q/Q[expected]),
               as.table=TRUE)

pdf(file.path(FOLDER,"Allplots","Qvalues-nFactors.pdf"),width=8,height=5)
print(xout)
dev.off()
