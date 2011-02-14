####################
## PMF execution and postprocessing program
## ~results_template.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


###_* load libraries
library(reshape)
library(fields)

###_* read matrix and wavenumbers
X <- read.table("matrix.dat")
xval <- scan("wavenumbers.txt",quiet=TRUE)

###_* functions

readF <- with(list(nc=length(xval)),
              function(x) matrix(scan(x,quiet=TRUE),ncol=nc,byrow=TRUE))

extractQval <- function(txt)
  as.numeric(sub(".*Final chi2=([ 0-9.]+).*","\\1",
                 grep("Final chi2",txt,value=TRUE)))

extractFPEAK <- function(txt)
  if(all(regexpr("fpeak",txt) < 1)) 0 else {
    as.numeric(sub('Parameter \\"fpeak\\" =(.*)',"\\1",
                   grep('Parameter \\"fpeak\\"',txt,value=TRUE)))
  }

lowerdiag <- function(M) sapply(1:ncol(M), function(i) M[i:nrow(M),i])

maybe <- function(x) if(length(x)==0) NA else x

f <- function(dnm) {
  F <- readF(file.path(dnm,"F_FACTOR.TXT"))
  txt <- readLines(file.path(dnm,"output.txt"))
  Qval <- maybe(extractQval(txt))
  fpeak <- maybe(extractFPEAK(txt))
  nFactors <- dim(F)[1]
  return(c(Qval=Qval,fpeak=fpeak,nFactors=nFactors))
}

g <- function(dnm) {
  F <- readF(file.path(dnm,"F_FACTOR.TXT"))
  G <- read.table(file.path(dnm,"G_FACTOR.TXT"))
  gcor <- local({x <- unlist(lowerdiag(cor(G)))
            median(x[x < 1.00])})
  txt <- readLines(file.path(dnm,"output.txt"))
  fpeak <- extractFPEAK(txt)
  nFactors <- dim(F)[1]
  return(c(gcor=gcor,fpeak=fpeak,nFactors=nFactors))
}

###_* apply functions

allf <- file.info(list.files("."))
dnames <- grep(".+\\_[0-9]{3}",rownames(allf)[allf$isdir],value=TRUE)
mymat <- as.data.frame(t(sapply(dnames,f)))

###_* plots

###_ . q-values

reshapen <- cast(mymat,fpeak~nFactors,value="Qval")
obj <- list(x=reshapen$fpeak,y=as.numeric(names(reshapen)[-1]),
            z = data.matrix(reshapen[,-1]))
pdf("Qvalues.pdf",paper="a4")
image.plot(obj,las=1,yaxt="n")
title(main="Q values",xlab="FPEAK",ylab="nFactors")
box()
par(new=TRUE)
plot(par("usr")[1:2],par("usr")[3:4],
     type="n",xlab="",ylab="",axes=FALSE,
     xaxs="i",yaxs="i")
axis(2,at=obj$y,las=1)
dev.off()

## interpolated version
#make.surface.grid( list(seq( min(obj$x),max(obj$x),,50),
#                        seq( min(obj$y),max(obj$y),,50)))-> loc
#interp.surface( obj, loc)-> look
#image.plot( as.surface( loc, look), las=1)
#box()
#title(xlab="FPEAK",ylab="nFactors")

###_ . correlations (G)

myGmat <- as.data.frame(t(sapply(dnames,g)))
reshapen <- cast(`names<-`(myGmat,sub("gcor","value",names(myGmat))),
                 fpeak~nFactors)
obj <- list(x=reshapen$fpeak,y=as.numeric(names(reshapen)[-1]),
            z = data.matrix(reshapen[,-1]))
pdf("Gmedian.pdf",paper="a4")
image.plot(obj,las=1,yaxt="n")
title(main="median correlations among all G",xlab="FPEAK",ylab="nFactors")
box()
par(new=TRUE)
plot(par("usr")[1:2],par("usr")[3:4],
     type="n",xlab="",ylab="",axes=FALSE,
     xaxs="i",yaxs="i")
axis(2,at=obj$y,las=1)
dev.off()
