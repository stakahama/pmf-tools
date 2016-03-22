####################
## PMF execution and postprocessing program
## ~fsmwa.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


## inputs:
input <- commandArgs()
pattern <- "--file=(.+)"
srcpath <- gsub('~+~'," ",dirname(sub(pattern,"\\1",input[grepl(pattern,input)])),fixed=TRUE)
source(file.path(srcpath,"functions/io.R"))

argv <- tail(input,-grep("--args",input,fixed=TRUE))
filename <- argv[1]

## contents
args <- readArgs(filename)
for(p in names(args))
  assign(p,args[[p]])

#########################################
# output: PercentInfo.pdf
#########################################
library(class)
library(fields)
X <- read.table(file.path(FOLDER,"matrix.dat"))
wn <- scan(file.path(FOLDER,"variables.txt"),quiet=TRUE)

## Full matrix
#pca.out <- prcomp(X,center=FALSE)
#svd.out <- svd(X)

#pca.factors <- t(pca.out$rotation)
#svd.factors <- t(svd.out$v)
#out <- knn(pca.factors,svd.factors,cl=1:nrow(pca.factors))
#svd.factors.sorted <- svd.factors[as.integer(as.character(out)),]
#svd.scores <- svd.out$d[as.integer(as.character(out))]

## Start FSWFA
n <- 20 ## window width in units of wavenumber indices
ind <- local({
  a <- 1:(ncol(X)-n)
  b <- a+n
  cbind(start=a,end=b)
})
fsw.scores <- apply(ind,1,function(.ind,.X)
                    rev(sort(svd(.X[,.ind[1]:.ind[2]])$d)),
                    .X=X)

W <- `rownames<-`((Eigenvalue=t(fsw.scores)),wn[ind[,1]+n/2])
Winfr <- t(apply(W,1,function(x) cumsum(x)/sum(x)))
r <- rbind(c(2380,2280),c(1550,1000))

library(fields)
dir.create(file.path(FOLDER,"Allplots"))
pdf(file.path(FOLDER,"Allplots","PercentInfo.pdf"))
local({
  X <- `colnames<-`(Winfr,1:ncol(Winfr))
  X <- X[order(as.numeric(rownames(X))),]

  xlim <- rev(range(as.numeric(rownames(X))))
  xlim <- c(ceiling(max(xlim)/1e2)*1e2,floor(min(xlim)/1e2)*1e2)
  par(mfrow=c(2,1),mar=c(4,4,1.5,6)+0.1,xaxs="i",yaxs="i",
      mgp=c(2.4,0.8,0),las=1)
                                        # plot 1
  matplot(as.numeric(rownames(W)),log10(W),type="l",lty=1,
          xlim=xlim,ylim=c(-5,0),
          xlab=expression("Wavenumber"~(cm^-1)),
          ylab=expression(log[10]("eigenvalue")),
          main="Eigenvalues")
  apply(r,1,function(x) rect(x[1],par("usr")[3],x[2],par("usr")[4],
                             col="white",border="white"))
  box()
                                        # plot 2
  image(x=as.numeric(rownames(X)),
        y=as.numeric(colnames(X)),
        z=X*100,axes=FALSE,
        col=tim.colors(64),
        xlim=xlim,ylim=c(15,1)+.5*c(1,-1),
        xlab=expression("Wavenumber"~(cm^-1)),
        ylab="Number of Factors",
        main="Percent Information")
  image.plot(z=range(X*100,na.rm=TRUE),legend.only=TRUE,
             col=tim.colors(64))
  xrg <- par("usr")[1:2]
  yrg <- par("usr")[3:4]
  apply(r,1,function(x) rect(x[1],yrg[1],x[2],yrg[2],
                             col="white",border="white"))
  abline(h=1:20 + 0.5,col="white",lty=0.9)
  box()
  par(new=TRUE)
  plot(NA,NA,xlim=xrg,ylim=yrg,
       axes=FALSE,type="n",xlab="",ylab="")
  axis(1)
  axis(2,at=1:ncol(Winfr),lab=1:ncol(Winfr),las=1,cex.axis=0.7,
       mgp=c(3,0.7,0))
})
dev.off()
