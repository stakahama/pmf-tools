####################
## PMF execution and postprocessing program
## ~stxmPMF.r~
## $Rev: -1 $
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


ev <-
  "278.00,278.50,279.00,279.50,280.00,280.50,281.00,281.50,282.00,282.50,283.00,283.50,284.00,284.50,284.65,284.80,284.95,285.10,285.25,285.40,285.55,285.70,285.85,286.00,286.15,286.30,286.45,286.60,286.75,286.90,287.05,287.20,287.35,287.50,287.65,287.80,287.95,288.10,288.25,288.40,288.55,288.70,288.85,289.00,289.15,289.30,289.45,289.60,289.75,289.90,290.05,290.20,290.35,290.50,290.65,290.80,290.95,291.10,291.25,291.40,291.55,291.70,291.85,292.00,292.50,293.00,293.50,294.00,294.50,295.00,295.50,296.00,296.50,296.80,297.09,297.38,297.67,297.96,298.26,298.55,298.84,299.13,299.42,299.71,300.00,302.00,304.00,306.00,308.00,310.00,312.00,314.00,316.00,318.00,320.00"
ev <- scan(textConnection(ev),sep=",",quiet=TRUE); gc()

library(TeachingDemos)
library(class)
library(reshape)
library(fields)
interp2D <- function(x,y,z) {
  obj <- list(x=x,y=y,z=z)
  return(obj)
  loc <-
    make.surface.grid(lapply(list(x=x,y=y),function(x)
                             do.call(seq,as.list(c(range(x),
                                                   length=length(x)*2)))))
  look <- interp.surface( obj, loc)
  as.surface(loc,look)
}

normalize <- function(x) (x-min(x))/diff(range(x))
rgbindex <- function(n)
  `names<-`(do.call(expand.grid,rep(list(seq(0,1,,round(n^(1/3)))),3)),
            c("red","green","blue"))
makergb <- function(U) do.call(rgb,unname(U))
cube2rgb <- function(rgbcube) do.call(rgb,rgbcube)
discretizeRGBcol <- function(lst,rgbcube)
  (function(rgbmat)
   as.integer(knn1(rgbmat,do.call(cbind,lst),1:nrow(rgbmat)))
   )(do.call(cbind,rgbcube))
normG <- function(fromfolder,bycol=TRUE) {
  G <- as.matrix(read.table(file.path(fromfolder,"G_FACTOR.TXT")))
  Gnorm <-
    (if( bycol ) lapply(split(G,col(G)),normalize)
    else (function(x) split(x,col(x)))(`[<-`(G,,value=normalize(c(G)))))
  if(length(Gnorm) < 3) Gnorm <- c(Gnorm,list(rep(0,length(Gnorm[[1]]))))
  Gnorm
}
## normG <- function(fromfolder,bycol=TRUE) {
##   G <- as.matrix(read.table(file.path(fromfolder,"G_FACTOR.TXT")))
##   #G[] <- apply(G,2,normalize)
##   Gnorm <- (function(x) split(x,col(x)))(sweep(G,1,rowSums(G),'/'))
##   if(length(Gnorm) < 3) Gnorm <- c(Gnorm,list(rep(0,length(Gnorm[[1]]))))
##   Gnorm
## }
getF <- function(fromfolder) 
        F <- as.matrix(read.table(file.path(fromfolder,"F_FACTOR.TXT")))

library(RSQLite)
dbname <- "../clusterspec/particledb.db"
drv <- dbDriver("SQLite")

##
n <- 64
con <- dbConnect(drv,dbname=dbname)
rgbcube <- rgbindex(n)
fi <- list.files(".","[_0-9]+p"); suffix <- c("002","006")
#for( f in (function(x) split(x,1:nrow(x)))(expand.grid(suffix,fi)[2:1]) ) {
for( f1 in fi ) {
  tb <- sub("([_0-9]+)p","P\\1",f1)
  print(tb)
  pngname <- sprintf("%s.png",sub("([0-9]+)p","\\1_pmf",f1))
  if(file.exists(file.path("figures",pngname))) next
  png(file.path("figures",pngname),width=96*6,height=96*5)
  par(mfrow=c(2,2),mar=c(4,4,1.5,1.5),mgp=c(2,1,0))
  for( f2 in suffix ) {
    result <- tryCatch({
      folder <- file.path(f1,paste(f1,f2,sep=""))
      if(all(!file.info(list.files(folder,"FACTOR",full=TRUE))$size > 0))
        return(NULL)
      xybw <- dbGetQuery(con,sprintf('SELECT x,y,bw FROM %s',tb))        
      ##
      Gnorm <- normG(folder)
      colind <- discretizeRGBcol(Gnorm,rgbcube)
      rc <- recast(cbind(xybw[1:2],
                         g=replace(rep(0,nrow(xybw)),xybw$bw==1,colind)),
                   x~y~variable,id.var=1:2)
      ##
      op <- (function(xy){
        (function(lims){
          op <- do.call(squishplot,c(lapply(xy,function(x) range(x)),1))
          plot.new(); do.call(plot.window,c(lims,xaxs="i",yaxs="i"))
          sapply(1:2,axis)
          do.call(rect,as.list(c(par("usr")[c(1,3,2,4)],col="dark blue")))
          do.call(image,c(do.call(interp2D,c(xy,list(z=rc[,,1]))),
                          list(col=cube2rgb(rgbcube),zlim=c(1,n),
                               ann=FALSE,asp=1,add=TRUE,
                               xaxs="i",yaxs="i")))
          abline(v=axTicks(1),col=8)
          abline(h=axTicks(2),col=8)
          title(xlab=expression("x"~(mu*m)),ylab=expression("y"~(mu*m)))
          op
        })(lapply(xy,range))
      })(lapply(dimnames(rc)[1:2],as.numeric))
      ##
      par(op)
      wholespec <-
        colMeans(dbGetQuery(con,
                            sprintf('SELECT * FROM %s WHERE bw = 1',
                                    tb))[-(1:3)])
      F <- t(getF(folder))
      plot(ev,wholespec,col=grey(0.45),lwd=3,type="l",
           xlab="ev",ylab="optical density (arb. units)")
      par(new=TRUE)      
      matplot(ev,F,lwd=1,type="l",lty=1,axes=FALSE,ann=FALSE,
              col=c("red","green","blue")[1:ncol(F)])
      1
    },error=function(e) NULL)
  }
  dev.off()
  if(is.null(result)) file.remove(file.path("figures",pngname))
}
dbDisconnect(con)
## matplot(t(F),type="l",lty=1,col=c("red","green","blue"))
## plot(1:8,col=myrgb(8),pch=15)

## particle <-
##   dbGetQuery(con,sprintf('SELECT x,y FROM %s WHERE bw=1',tb))       
## with(particle,plot(x,y,col=do.call(rgb,unname(Gnorm)),
##                    pch=19,cex=2.5))
## with(particle,identify(x,y))
## image(matrix(sample(c(49,49,4),40,rep=TRUE),ncol=4),zlim=c(1,64),
##       col=cube2rgb(rgbcube))
## a <- discretizeRGBcol(lapply(Gnorm,`[`,c(71,98,813)),rgbcube)
## rgbcube[a,]

## plot(1:78,col=do.call(rgb,unname(as.data.frame(Gnorm)[which(with(particle, x > 7 & x < 9 & y > 5 & y < 7)),])),pch=15)
## plot(1:64,col=myrgb(64),pch=15,cex=2)



