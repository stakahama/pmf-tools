###_* user inputs

source("userinputs.r")

###_* functions

library(TeachingDemos)

toimage <- function(pixels,ev,mycol) {
   nfac <- ncol(ev)-1
   pixcol <- do.call(rgb,replace(rep(list(0),3),1:nfac,unname(ev[,1:nfac])))
   ##
   im <- reshape(data.frame(pixels,z=as.integer(factor(pixcol))),
                 timevar="y",idvar="x",direction="wide")
   y <- as.numeric(substring(names(im)[-1],3))
   j <- order(y)
   ##
   list(x = as.numeric(im$x),
        y = y[j],
        z = as.matrix(im[,-1])[,j],
        col = pixcol)
}

###_* inputs

energies <- scan(file.path(FOLDER,"variables.txt"),quiet=TRUE)
pixels <- read.table(file.path(FOLDER,"samples.txt"),
                     colClasses="character",
                     col.names=c("x","y"),sep=",")
simgrid <- read.delim(file.path(FOLDER,"simgrid.txt"),row.names=1)

###_* plots

specpath <- file.path(FOLDER,"Allplots","spectra")
dir.create(specpath)
for( runpath in list.files(FOLDER,basename(FOLDER),full=TRUE) ) {
  ff <- as.matrix(read.table(file.path(runpath,"F_FACTOR.txt")))
  if( nrow(ff) > 3 ) next()
  ev <- read.delim(file.path(runpath,"ExplainedVariation.txt"))
  im <- toimage(pixels,ev)
  png(file.path(specpath,sprintf("%s.png",basename(runpath))),
      width=8*96,height=5*96)
  par(mfrow=c(1,2),mar=c(3.5,3.5,.5,.5),mgp=c(2.2,.5,0),oma=c(0,0,1.5,0))
  matplot(energies,t(ff),type="l",lty=1,
          col=c("red","green","blue")[1:nrow(ff)],
          xlab="Energy (eV)",
          ylab="Absorbance (arb. units)")
  squishplot(range(im$x),range(im$y),asp=1,newplot=FALSE)
  with(im,image(x,y,z,col=col,
                xlab=expression("x"~(mu*m)),
                ylab=expression("y"~(mu*m))))
  mtext(paste(sprintf("%s=%s",names(simgrid),simgrid[basename(runpath),]),
              collapse=", "),3,outer=TRUE)
  dev.off()
}

