####################
## PMF execution and postprocessing program
## ~errormatrix.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################

library(fields)
library(Hmisc)
library(TeachingDemos)

source("userinputs.r")
source("functions/classify.r")
source("functions/imageplots.r")

Arg <- tail(commandArgs(),1)
runno <- as.integer(Arg)
if( is.na(runno) )
  stop("enter run number as integer")

solution <- runnum(runno)

mz <- scan(file.path(FOLDER,"variables.txt"),quiet=TRUE)
samples <- scan(file.path(FOLDER,"samples.txt"),"",sep="\n",quiet=TRUE)
X <- as.matrix(read.table(file.path(FOLDER,"matrix.dat")))
S <- as.matrix(read.table(file.path(FOLDER,"std_dev.dat")))
F <- as.matrix(read.table(file.path(solution,"F_FACTOR.TXT")))
G <- as.matrix(read.table(file.path(solution,"G_FACTOR.TXT")))

rownames(X) <- rownames(G) <- samples
rownames(F) <- colnames(G) <- sprintf("Factor%d",1:nrow(F))
E <- X - G%*%F

##{{{
## E <- X - G%*%F
## Cx <- 1/(nrow(E)-1)*t(E)%*%E
## Sx <- Cx/with(list(s=apply(E,2,sd)),outer(s,s,`*`))
## Ecov <- cov(E)
## Ecor <- cor(E)
##

## par(mfrow=c(2,2),mar=c(4,4,1,1))
## with(list(lim=range(Ecov,Cx)),{
##   squishplot(lim,lim,asp=1,new=FALSE)
##   plot(Ecov,Cx,asp=1,pch=19,cex=.5,xlim=lim,ylim=lim)
##   abline(0,1)  
## })
## with(list(lim=c(-1,1)),{
##   squishplot(lim,lim,asp=1,new=FALSE)
##   plot(Ecor,Sx,asp=1,xlim=lim,ylim=lim,pch=19,cex=.5,)
##   abline(0,1)  
## })
## with(list(lim=range(mz)),{
##   squishplot(lim,lim,asp=1,new=FALSE)
##   image(mz,mz,Sx,zlim=c(-1,1),asp=1,zlim=c(-1,1))
##   lapply(1:2,axis)
##   squishplot(lim,lim,asp=1,new=FALSE)  
##   image(mz,mz,Ecor,zlim=c(-1,1),asp=1,zlim=c(-1,1))
##   lapply(1:2,axis)  
## })
##}}}

filename <-
  with(list(s=sprintf("%s_%03d-PMF-Residuals-%%s%%s",basename(FOLDER),runno)),
       function(x,ext=".png") sprintf(s,x,ext))

png(file.path(FOLDER,"Allplots",filename("E")),width=96*6)
par(mar=c(4,5,1.5,6),mgp=c(2.2,.5,0))
image(mz,rev(1:length(samples)),t(E),col=tim.colors(64),
      xaxs="i",yaxs="i")
axis(1)
with(list(s=samples),{
  ix <- seq(1,length(s),floor(length(s)/20))  
  axis(2,at=ix,s[ix],cex.axis=.5,las=1)
})
box()
title(xlab="m/z")
mtext("E",3)
image.plot(z=t(E),col=tim.colors(64),legend.only=TRUE)
dev.off()

png(file.path(FOLDER,"Allplots",filename("ES")),width=96*6)
par(mar=c(4,5,1.5,6),mgp=c(2.2,.5,0))
image(mz,rev(1:length(samples)),t(E/S),col=tim.colors(64),
      xaxs="i",yaxs="i")
axis(1)
with(list(s=samples),{
  ix <- seq(1,length(s),floor(length(s)/20))
  axis(2,at=ix,s[ix],cex.axis=.5,las=1)
})
box()
title(xlab="m/z")
mtext("E/S",3)
image.plot(z=t(E/S),col=tim.colors(64),legend.only=TRUE)
dev.off()

Ecor <- cor(E)
n <- 4

png(file.path(FOLDER,"Allplots",filename("mz")),width=96*5.5)
par(mar=c(4,4,1.5,6),mgp=c(2.2,.5,0))
squishplot(range(mz),range(mz),asp=1)
image(mz,mz,t(Ecor),col=tim.colors(n),
      xaxs="i",yaxs="i",asp=1,zlim=c(-1,1))
axis(1)
axis(2)
box()
title(xlab="m/z",ylab="m/z")
mtext(expression("normalized"~ E^T*E),3,line=.5)
image.plot(zlim=c(-1,1),col=tim.colors(n),
           legend.only=TRUE)
dev.off()

plotbox <- function(mat,mz,ev=FALSE) {
  plot.new()
  if(ev) {
    plot.window(range(mz),c(0,100))
  } else {
    plot.window(range(mz),
                quantile(mat,.995,na.rm=TRUE)*c(-1,1))
  }
  abline(h=0,col=8)
  boxplot(split(mat,col(mat)), xaxt="n",lty=1,
          pars=list(outpch=".",outcol="blue"),
          border="blue",add=TRUE,at=mz,lend=3)
  axis(1,pretty(mz))
}

pdf(file.path(FOLDER,"Allplots",filename("boxplots",".pdf")),
    width=8,height=6)
par(mfrow=c(3,1))
par(mar=c(1,1,.5,.5),oma=c(3,3,1,1),mgp=c(2.2,.5,0))    
plotbox(E,mz)
title(ylab="E",xpd=NA)
plotbox(E/S,mz)
title(ylab="E/S",xpd=NA)
plotbox(abs(G%*%F)/(abs(G%*%F)+abs(E))*1e2,mz,TRUE)
axis(1,pretty(mz))
title(ylab="Explained Variation (%)",xpd=NA)
mtext("m/z",1,outer=TRUE,line=1)
dev.off()

EVi <- rowSums(abs(G%*%F)/S)/
  rowSums((abs(G%*%F)+abs(E))/S)*1e2
Org <- rowSums(X[,na.omit(match(10:125,mz))])

pdf(file.path(FOLDER,"Allplots",filename("EV-Org",".pdf")),
    width=8,height=5)
plot(Org,EVi,
     xlim=c(0,max(Org,na.rm=TRUE)*1.2),
     ylim=c(0,100),
     col=4,cex=.5,
     xlab=expression("Org"~(mu*g/m^3)),
     ylab="Explained Variation (%)")
dev.off()
