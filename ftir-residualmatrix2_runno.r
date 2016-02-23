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


input <- commandArgs()
pattern <- "--file=(.+)"
srcpath <- gsub('~+~'," ",dirname(sub(pattern,"\\1",input[grepl(pattern,input)])),fixed=TRUE)
source(file.path(srcpath,"functions/io.R"))

argv <- tail(input,-grep("--args",input,fixed=TRUE))
filename <- argv[1]
runno <- as.integer(argv[2])

if( is.na(runno) )
  stop("enter run number as integer")

args <- read.args(filename)
for(p in names(args))
  assign(p,args[[p]])

source("functions/classify.r")
source("functions/imageplots.r")

extendimage.symmetric <- function(x,z) {
  newx <- seq(head(x,1),tail(x,1),by=median(diff(x)))
  newz <- matrix(NA,length(newx),length(newx))
  ix <- sapply(x,function(.x,y) which.min(abs(.x-y)),newx)
  newz[ix,ix] <- z
  list(xy=newx,z=newz)
}

m <- extendimage.symmetric(var,Ecor)
with(m,filled.contour(xy,xy,z))
with(m,image(t(z)))

solution <- runnum(runno)

var <- scan(file.path(FOLDER,"variables.txt"),quiet=TRUE)
samples <- scan(file.path(FOLDER,"samples.txt"),"",sep="\n",quiet=TRUE)
X <- as.matrix(read.table(file.path(FOLDER,"matrix.dat")))
S <- as.matrix(read.table(file.path(FOLDER,"std_dev.dat")))
F <- as.matrix(read.table(file.path(solution,"F_FACTOR.TXT")))
G <- as.matrix(read.table(file.path(solution,"G_FACTOR.TXT")))

rownames(X) <- rownames(G) <- samples
rownames(F) <- colnames(G) <- sprintf("Factor%d",1:nrow(F))
E <- X - G%*%F

filename <-
  with(list(s=sprintf("%s_%03d-PMF-Residuals-%%s%%s",basename(FOLDER),runno)),
       function(x,ext=".png") sprintf(s,x,ext))

Ecor <- cor(E)
n <- 4
lim <- rev(range(var))

png(file.path(FOLDER,"Allplots",filename("Wavenumber")),width=96*5.5)
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,6),mgp=c(2.2,.5,0))
squishplot(lim,lim,asp=1)
plot.window(lim,lim,xaxs="i",yaxs="i")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="black")
with(extendimage.symmetric(var,Ecor),
     image(rev(xy),rev(xy),t(z),col=tim.colors(n),
           zlim=c(-1,1),add=TRUE))
axis(1)
axis(2)
box()
title(xlab="Wavenumber",ylab="Wavenumber")
mtext(expression("normalized"~ E^T*E),3,line=.5)
image.plot(zlim=c(-1,1),col=tim.colors(n),
           legend.only=TRUE)
dev.off()

plotbox <- function(mat,var,ev=FALSE,xlim=c(4000,1200)) {
  plot.new()
  if(ev) {
    plot.window(xlim,c(0,100))
  } else {
    plot.window(xlim,quantile(abs(mat),.999)*c(-1,1))
  }
  abline(h=0,col=8)
  qy <- apply(mat,2,quantile,c(.05,.25,.5,.75,.95))
  with(structure(split(qy,row(qy)),names=rownames(qy)), {
    polygon(c(var,rev(var)),c(`5%`,rev(`95%`)),border=NA,
            col=rgb(t(col2rgb(grey(.7)))/255,alpha=.5))
    polygon(c(var,rev(var)),c(`25%`,rev(`75%`)),border=NA,
            col=rgb(t(col2rgb(grey(.5)))/255,alpha=.5))
    lines(var,`50%`,col=4)
  })
  axis(1,pretty(var))
  axis(2)
  box()  
}

## smoothScatter(rep(var,nrow(E)),c(t(E)),xlim=c(4000,1200),
##               xlab="Wavenumber",ylab="E")
## smoothScatter(rep(var,nrow(E)),c(t(E/S)),xlim=c(4000,1200),
##               xlab="Wavenumber",ylab="E/S")
## smoothScatter(rep(var,nrow(E)),1e2*c(t(abs(G%*%F)/(abs(G%*%F)+abs(E)))),
##               xlim=c(4000,1200),
##               xlab="Wavenumber",ylab="Explained Variatoin (%)")

pdf(file.path(FOLDER,"Allplots",filename("boxplots",".pdf")),
    width=8,height=6)
par(mfrow=c(3,1))
par(mar=c(1,1,.5,.5),oma=c(3,3,1,1),mgp=c(2.2,.5,0))    
plotbox(E,var)
title(ylab="E",xpd=NA)
plotbox(E/S,var)
title(ylab="E/S",xpd=NA)
plotbox(abs(G%*%F)/(abs(G%*%F)+abs(E))*1e2,var,TRUE)
title(ylab="Explained Variation (%)",xpd=NA)
mtext("Wavenumber",1,outer=TRUE,line=1)
dev.off()

EVi <- rowSums(abs(G%*%F)/S)/
  rowSums((abs(G%*%F)+abs(E))/S)*1e2
area <- apply(X,1,function(x,d) sum(0.5*(head(x,-1)+tail(x,-1))*d),abs(diff(var)))

pdf(file.path(FOLDER,"Allplots",filename("EV-Org",".pdf")),
    width=8,height=5)
plot(area,EVi,
     xlim=c(0,max(area,na.rm=TRUE)*1.2),
     ylim=c(0,100),
     col=4,cex=.5,
     xlab=expression("Integrated area"~(abs)),
     ylab="Explained Variation (%)")
dev.off()
