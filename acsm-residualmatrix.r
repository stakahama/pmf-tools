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
  with(list(s=sprintf("PMF-Residuals-%%s_%s_%03d.png",basename(FOLDER),runno)),
       function(x) sprintf(s,x))

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
mtext("E",3)
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
