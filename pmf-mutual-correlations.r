###_* load libraries and functions

library(lattice)
library(latticeExtra)

read.F <- function(x) {
  ff <- as.matrix(read.table(x))
  ctabl <- cor(t(ff))
  clabl <- outer(sprintf("%02d",1:nrow(ctabl)),
                 sprintf("%02d",1:ncol(ctabl)),
                 paste,sep="-")
  structure(ctabl[upper.tri(ctabl)],
            names=clabl[upper.tri(clabl)])
}
read.G <- function(x) {
  gg <- as.matrix(read.table(x))
  ctabl <- cor(gg)
  clabl <- outer(sprintf("%02d",1:nrow(ctabl)),
                 sprintf("%02d",1:ncol(ctabl)),
                 paste,sep="-")
  structure(ctabl[upper.tri(ctabl)],
            names=clabl[upper.tri(clabl)])
}

###_* inputs
source("userinputs.r")

simgrid <- read.delim(file.path(FOLDER,"simgrid.txt"),row.names=1)
for( x in 1:ncol(simgrid) ) simgrid[,x] <- factor(simgrid[,x])
runs <- list.files(FOLDER,basename(FOLDER),full=TRUE)

###_* calculation of correlations
allcor <- do.call(rbind,lapply(runs,function(x) {
  tryCatch({
    gcor <- read.G(file.path(x,"G_FACTOR.TXT"))
    data.frame(run=basename(x),
               pairs=names(gcor),
               G=gcor,
               F=read.F(file.path(x,"F_FACTOR.TXT"))[names(gcor)])
  },error=function(e) NULL)
}))
mg <- merge(allcor,simgrid,,1,0,all=TRUE)
write.table(mg,file.path(FOLDER,"Allplots","mutual-correlations.txt"),
            sep="\t",row.names=FALSE,quote=FALSE)

###_* plot
out <- xyplot(G~F | FPEAK*nFactors, groups=Seed,
              data=mg,
              panel=function(...)
              { panel.abline(v=0,lty=3,col=8)
                panel.abline(h=0,lty=3,col=8)
                panel.xyplot(...)
              },
              scales=list(at=seq(-1,1,.5)),
              xlim=c(-1,1),ylim=c(-1,1),
              auto=list(space="right",title="Seed",cex.title=.9),
              par.settings=list(superpose.symbol=list(pch=c(1,3,4),cex=.6)),
              as.table=TRUE,
              xlab="F correlation",ylab="G correlation")
pdf(file.path(FOLDER,"Allplots","mutual-correlations.pdf"),width=12,height=8)
print(useOuterStrips(out))
dev.off()

###_* example identification

## mg <- read.delim(file.path(FOLDER,"Allplots","mutual-correlations.txt"))

## with(subset(mg,FPEAK==-1.2 & Seed==100 & nFactors==3),{
##   plot(F,G,xlim=c(-1,1),ylim=c(-1,1),pch=3,col=4,lwd=2,ann=FALSE)
##   abline(h=0,lty=3,col=8);abline(v=0,lty=3,col=8)
##   text(F,G,adj=c(0,1),pairs,xpd=NA)
##   title(xlab="F correlation",ylab="G correlation",
##         main=paste(sprintf("%s = %s",c("nFactors","FPEAK","Seed"),
##           c(nFactors[1],sprintf("%.1f",FPEAK[1]),Seed[1])),collapse=", "))
##   ##identify(F,G,pairs)
## })
