####################
## PMF execution and postprocessing program
## ~converg.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


allfiles <- list.files(".",patt="output\\.txt",rec=TRUE,full=TRUE)

mtime <-
  `names<-`(sapply(allfiles,function(x) {
    file.info(x)[["mtime"]] > ISOdate(2007,12,31,13,20,tz="PST")
  }),substring(dirname(allfiles),3))
  
conv <- `names<-`(sapply(allfiles,function(x)
                     all(regexpr("NO CONVERGENCE",readLines(x))<0)),
              substring(dirname(allfiles),3))

#as.matrix(conv)
abc <- data.frame(mtime=as.matrix(mtime),conv=as.matrix(conv))
#subset(abc,mtime)

simgrid <- read.table("simgrid.txt")
mydf <- local({x <- merge(simgrid,abc[,"conv",drop=FALSE],by=0,all=TRUE);
               `rownames<-`(x[,-1],x[,1])})

library(reshape)
library(fields)
reshapen <- cast(`names<-`(mydf,sub("conv","value",names(mydf))),
                 FPEAK~nFactors)
obj <- list(x=reshapen$FPEAK,y=as.numeric(names(reshapen)[-1]),
      z = data.matrix(reshapen[,-1]))
pdf("convergence.pdf")
image.plot(obj,col=c("red","midnightblue"),yaxt="n",
           axis.args=list(at=c(0,1),lab=c("False","True")))
title(main="Convergence",xlab="FPEAK",ylab="nFactors")
par(new=TRUE)
plot.window(xlim=range(obj$x),ylim=range(obj$y))
axis(2,at=obj$y,las=1)
box()
dev.off()
