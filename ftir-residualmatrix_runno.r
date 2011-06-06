####################
## PMF execution and postprocessing program
## ~errormatrix.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################

source("userinputs.r")
source("functions/classify.r")
source("functions/imageplots.r")

Arg <- tail(commandArgs(),1)
runno <- as.integer(Arg)
if( is.na(runno) )
  stop("enter run number as integer")

solution <- runnum(runno)

wn <- scan(file.path(FOLDER,"variables.txt"),quiet=TRUE)
samples <- scan(file.path(FOLDER,"samples.txt"),"",quiet=TRUE)
X <- as.matrix(read.table(file.path(FOLDER,"matrix.dat")))
S <- as.matrix(read.table(file.path(FOLDER,"std_dev.dat")))
F <- as.matrix(read.table(file.path(solution,"F_FACTOR.TXT")))
G <- as.matrix(read.table(file.path(solution,"G_FACTOR.TXT")))

rownames(X) <- rownames(G) <- samples
rownames(F) <- colnames(G) <- sprintf("Factor%d",1:nrow(F))

E <- X - G%*%F
Ecov <- E%*%t(E)
Ecov <- cor(E)

library(fields)
##{{{
## image.plot(t(Ecov[1,,drop=FALSE]))
## image.plot(Ecov[1:5,],axes=FALSE)
## axis(1,axTicks(1),lab=approx(1:ncol(E),wn,axTicks(1))$y)
##}}}
library(Hmisc)
##{{{
## par(mar=c(4,4,2,6))
## xlabs <- replace(seq(4000,1500,-500),1,3800)
## xticks <- approxExtrap(wn,1:ncol(E),xlabs)$y
## image(1:ncol(E),1:nrow(E),t(E),col=tim.colors(64),axes=FALSE,
##       xlim=range(xticks),ann=FALSE)
## axis(1,tail(xticks,-1),lab=tail(xlabs,-1),mgp=c(2.2,.7,0))
## axis(2,1:nrow(E),lab=rownames(E),las=1,cex.axis=.5,tck=-.01,mgp=c(2.2,.5,0))
## box()
## title(xlab=expression(Wavenumber ~ (cm^-1)),mgp=c(2.2,.7,0))
## title(ylab="Sample",mgp=c(2.5,.7,0))
## image.plot(t(E),col=tim.colors(64),legend.only=TRUE)
##}}}

extwn <- seq(max(wn),min(wn),median(diff(wn)))
asint <- function(x) as.integer(round(x*100))
extE <- do.call(cbind,
                Map(function(w,x,e,n)
                    if(w %in% x) e[,x %in% w] else rep(NA,n),
                    asint(extwn), MoreArgs=list(x=asint(wn), e=E, n=nrow(E))))

extEs <- do.call(cbind,
                Map(function(w,x,e,n)
                    if(w %in% x) e[,x %in% w] else rep(NA,n),
                    asint(extwn), MoreArgs=list(x=asint(wn), e=E/S,
                                    n=nrow(E))))


library(TeachingDemos)

repl <- function(x) {
  i <- is.na(x)
  x[i] <- 1
  x[!i] <- NA
  x
}

### ---

filename <-
  with(list(s=sprintf("PMF-Residuals-%%s_%s_%03d.pdf",basename(FOLDER),runno)),
       function(x) sprintf(s,x))

pdf(file.path(FOLDER,"Allplots",filename("E")))
##
xlabs <- replace(seq(4000,1500,-500),1,3800)
xticks <- approxExtrap(extwn,1:ncol(extE),xlabs)$y
ylim <- c(1,nrow(extE))
par(mar=c(4,4,2,6))
##{{{---black strips---
## image(1:ncol(extE),1:nrow(extE),repl(t(extE)),col="black",
##       xlim=range(xticks),ylim=ylim,axes=FALSE,ann=FALSE)
##}}}
plot.new()
plot.window(range(xticks),ylim,xaxs="i",yaxs="i")
do.call(rect,c(as.list(par("usr")[c(1,3,2,4)]),col="black"))
image(1:ncol(extE),1:nrow(extE),t(extE),col=tim.colors(64),add=TRUE)
axis(1,tail(xticks,-1),lab=tail(xlabs,-1),mgp=c(2.2,.7,0),xpd=TRUE)
axis(2,1:nrow(E),lab=rownames(E),las=1,cex.axis=.5,tck=-.01,mgp=c(2.2,.5,0))
box()
title(xlab=expression(Wavenumber ~ (cm^-1)),mgp=c(2.2,.7,0))
title(ylab="Sample",mgp=c(2.5,.7,0))
title(main="Residual Matrix, E")
image.plot(t(extE),col=tim.colors(64),legend.only=TRUE)
##
dev.off()

pdf(file.path(FOLDER,"Allplots",filename("ES")))
local({
  extE <- extEs
  xlabs <- replace(seq(4000,1500,-500),1,3800)
  xticks <- approxExtrap(extwn,1:ncol(extE),xlabs)$y
  ylim <- c(1,nrow(extE))
  par(mar=c(4,4,2,6))
  plot.new()
  plot.window(range(xticks),ylim,xaxs="i",yaxs="i")
  do.call(rect,c(as.list(par("usr")[c(1,3,2,4)]),col="black"))
  image(1:ncol(extE),1:nrow(extE),t(extE),col=tim.colors(64),add=TRUE)
  axis(1,tail(xticks,-1),lab=tail(xlabs,-1),mgp=c(2.2,.7,0),xpd=TRUE)
  axis(2,1:nrow(E),lab=rownames(E),las=1,cex.axis=.5,tck=-.01,mgp=c(2.2,.5,0))
  box()
  title(xlab=expression(Wavenumber ~ (cm^-1)),mgp=c(2.2,.7,0))
  title(ylab="Sample",mgp=c(2.5,.7,0))
  title(main="Scaled Residual Matrix, E/S")
  image.plot(t(extE),col=tim.colors(64),legend.only=TRUE)
})
dev.off()

pdf(file.path(FOLDER,"Allplots",filename("Wavenumbers")))
local({
  par(mar=c(4,4,2,6))  
  extE <- cor(extE)
  xlabs <- replace(seq(4000,1500,-500),1,3800)
  xticks <- approxExtrap(extwn,1:ncol(extE),xlabs)$y
  xlim <- c(1,ncol(extE))
  squishplot(xlim=xlim,ylim=xlim,asp=1,new=FALSE)  
  plot.new()
  plot.window(xlim=xlim,ylim=xlim,asp=1,xaxs="i",yaxs="i")
  do.call(rect,c(as.list(par("usr")[c(1,3,2,4)]),col="black"))
  image(1:ncol(extE),1:nrow(extE),t(extE),col=tim.colors(64),zlim=c(-1,1),
        add=TRUE)
  axis(1,tail(xticks,-1),lab=tail(xlabs,-1),mgp=c(2.2,.7,0),xpd=TRUE)
  axis(2,tail(xticks,-1),lab=tail(xlabs,-1),mgp=c(2.2,.7,0),xpd=TRUE)
  box()
  title(xlab=expression(Wavenumber ~ (cm^-1)),mgp=c(2.2,.7,0))
  title(ylab=expression(Wavenumber ~ (cm^-1)),mgp=c(2.2,.7,0))
  title(main=expression("normalized"~E^T*E))
  image.plot(t(extE),col=tim.colors(64),zlim=c(-1,1),legend.only=TRUE)
})
dev.off()

pdf(file.path(FOLDER,"Allplots",filename("Samples")))
local({
  par(mar=c(4,4,2,6))  
  E <- cor(t(E))
  xticks <- 1:length(samples)
  squishplot(xlim=range(xticks),ylim=range(xticks),asp=1,new=FALSE)
  image(1:ncol(E),1:nrow(E),t(E),col=tim.colors(64),zlim=c(-1,1),axes=FALSE,
        xlim=range(xticks),ylim=range(xticks),ann=FALSE,asp=1,
        xaxs="i",yaxs="i")
  axis(1,xticks,lab=FALSE,mgp=c(2.2,.5,0),las=1,tck=-.01)
  text(xticks,par("usr")[3]-par("cxy")[2]*1.1,samples,cex=.5,xpd=NA,srt=90)
  axis(2,xticks,lab=samples,mgp=c(2.2,.5,0),las=1,tck=-.01,cex.axis=.5)
  box()
  title(xlab="Samples",mgp=c(2.5,.7,0))
  title(ylab="Samples",mgp=c(2.5,.7,0))
  title(main=expression("normalized"~E*E^T))  
  image.plot(t(E),col=tim.colors(64),zlim=c(-1,1),legend.only=TRUE)
})
dev.off()

