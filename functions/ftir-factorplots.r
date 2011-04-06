####################
## PMF execution and postprocessing program
## ~factorplots.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


###_* input

allf <- file.info(list.files("."))
dnames <- grep(".+\\_[0-9]{3}",rownames(allf)[allf$isdir],value=TRUE)
simgrid <- read.table("simgrid.txt")
xval <- scan('variables.txt',what=0,quiet=TRUE)
readF <- with(list(nc=length(xval)),
              function(x) matrix(scan(x,quiet=TRUE),ncol=nc,byrow=TRUE))

###_* functions

###_ . make one plot

## doPlot <- function(i,F,xval) {
##   ## for STXM
##   plot(xval,F[i,],type="l",lwd=2,col="dark blue")
## }
doPlot <- local({
  wn <- xval
  wndif <- round(diff(wn),2)
  winds <- local({j <- which(abs(wndif) > abs(median(wndif)))
                  t(matrix(sort(c(max(wn),wn[j],wn[j+1],min(wn))),nrow=2))})
  specs <- list(alcohol=3300,
                alkane=c(2925,2882,2852),
                alkene=2980,
                aromatic=3050,
                carbonyl=1720,
                organoS=876)
  function(i,F,wn) {
    ## for FTIR
    plot.new()
    plot.window(xlim=c(3800,1200),ylim=c(0,quantile(F,0.98)))
    title(xlab="",ylab="")
    abline(v=unlist(specs),col="orange")
    apply(winds,1,function(x)
          {ind <- (wn > x[1] & wn < x[2])
           lines(wn[ind],F[i,ind],col=4,lwd=2)})
    axis(1,cex.axis=0.8,mgp=c(1,0.5,0),xpd=NA)
    invisible(lapply(2:4,axis,label=FALSE))
    title(ylab="Absorbance",mgp=c(1.5,0,0))
    box()
  }
})

###_ . make lots of plots

readPlot <- function(dnm,xval,simgrid) {
  peq <- function(x) paste(x,collapse=" = ")
  lab <- unlist(rbind(colnames(simgrid),simgrid[dnm,]))  
  F <- readF(file.path(dnm,"F_FACTOR.TXT"))
  pdf(file.path("Allplots",paste(dnm,".pdf",sep="")),
      paper="a4",width=6.5,height=9)
  if( nrow(F) > 0 ) {
    par(mfrow=n2mfrow(nrow(F)),mar=c(2,1,1,1),oma=c(0,0,3,0),pty="s")
    sapply(1:nrow(F),doPlot,F,xval)
  } else {
    par(mfrow=c(4,4),mar=c(2,1,1,1),oma=c(0,0,3,0),pty="s")
    sapply(1:16,function(x) plot.new())
  }
  mtext(paste(dnm,peq(lab[1:2]),peq(lab[3:4]),
              peq(lab[5:6]),sep=", "),3,outer=TRUE)
  dev.off()
}

###_* apply

invisible(sapply(dnames,readPlot,xval,simgrid))


