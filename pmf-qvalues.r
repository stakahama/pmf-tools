####################
## PMF execution and postprocessing program
## ~qvalues.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


###_* command-line arguments
Arg <- tail(commandArgs(),1)
if( Arg=="--args" || Arg=="1" ) {
  CHI2 <- "final"
  OUTFILE1 <- "Qvalues-FPEAK.pdf"
  OUTFILE2 <- "Qvalues-nFactors.pdf"
} else {
  CHI2 <- "true"
  OUTFILE1 <- "Qvalues_true-FPEAK.pdf"
  OUTFILE2 <- "Qvalues_true-nFactors.pdf"
}

###_* read simgrid

###_ . inputs
source("userinputs.r")
simgrid <- read.delim(file.path(FOLDER,"simgrid.txt"),row.names=1,as.is=TRUE)
simgrid <- simgrid[with(simgrid,order(nFactors,FPEAK,Seed)),]
dir.create(file.path(FOLDER,"Allplots"))

###_* read dim(X), Q-values

###_ . X
Xmat <- read.table(file.path(FOLDER,"MATRIX.DAT"))
eqm <- function(p) prod(dim(Xmat)) - p*(sum(dim(Xmat)))

###_ . Q-values
outputs <- list.files(FOLDER,"output.txt",recursive=TRUE,full=TRUE)
maybe <- function(x) if(length(x)==0) NA else x
if( CHI2=="final" ) {
  qval <- function(x, patt=".*Final chi2[=][ ]*([0-9.]+).*")
    maybe(as.numeric(sub(patt,"\\1",grep(patt,readLines(x),value=TRUE))))
} else {
  qval <- function(x, patt="The true final chi2 value.+[=][ ]*([0-9.]+).*")
    maybe(as.numeric(sub(patt,"\\1",grep(patt,readLines(x),value=TRUE))))
}
vals <- sapply(outputs,qval)
names(vals) <- basename(dirname(names(vals)))

###_* plot
library(lattice)
simgrid$`Q-ratio` <- vals[rownames(simgrid)]/eqm(simgrid$nFactors)
for( var in c("nFactors","Seed") )
  simgrid[,var] <- factor(simgrid[,var])

panel <- function(...) {
  panel.abline(h=0,col=8)
  panel.xyplot(...)
}
qvalplot <- function(yaxis.relation="sliced")  {
  simgrid <- simgrid[order(simgrid$FPEAK),]
  xyplot(`Q-ratio`~FPEAK | nFactors, groups=Seed,data=simgrid,
         panel=panel,scales=list(y=list(relation=yaxis.relation)),
         type="o",auto=list(space="right",title="Seed",cex.title=1),
         ylab=expression(Q/Q[expected]),
         as.table=TRUE)
}
xout <- tryCatch(qvalplot("sliced"),error=function(e)
                 qvalplot("same"))

pdf(file.path(FOLDER,"Allplots",OUTFILE1),width=8,height=5)
## print(update(xout,scales=list(y=list(relation="same"))))
print(xout)
dev.off()

panel <- function(x,y,...,pch,groups,subscripts) {
  panel.points(x,y,pch=pch[subscripts],
               col=trellis.par.get("superpose.symbol")$col[groups[subscripts]])
  tryCatch(panel.lines(smooth.spline(x,y)),error=function(e) panel.lmline(x,y))
}

simgrid$FPEAK <- factor(round(simgrid$FPEAK,2))
mykey <- with(simgrid,{
  p <- trellis.par.get("superpose.symbol")
  dfr <- rbind(data.frame(lab=paste("Seed =",levels(Seed)),col=p$col[1:nlevels(Seed)],pch=1),
               data.frame(lab=paste("FPEAK =",levels(FPEAK)),col=1,pch=1:nlevels(FPEAK)))
  list(space="right",text=dfr["lab"],points=dfr[c("col","pch")])
})

mykey <- with(simgrid,{
  p <- trellis.par.get("superpose.symbol")
  lab <- do.call(paste,c(sep=",",rev(expand.grid(levels(Seed),levels(FPEAK)))))
  pp <- rev(expand.grid(col=p$col[1:nlevels(Seed)],pch=1:nlevels(FPEAK),
                    stringsAsFactors=FALSE))
  list(space="right",title="FPEAK,Seed",cex.title=.9,
       text=list(lab=lab),points=pp)
})

qvalnplot <- function(relation)
  xyplot(`Q-ratio`~nFactors, groups=Seed,data=simgrid,
         panel=panel,pch=simgrid$FPEAK,
         scales=list(y=list(relation=relation)),
         key=mykey,
         ylab=expression(Q/Q[expected]),
         as.table=TRUE)

simgrid$nFactors <- as.integer(as.character(simgrid$nFactors))
xout <- tryCatch(qvalnplot("sliced"),error=function(e)
                 qvalnplot("same"))


pdf(file.path(FOLDER,"Allplots",OUTFILE2),width=8,height=5)
print(xout)
dev.off()

