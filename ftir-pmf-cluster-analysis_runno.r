###_1. load libraries
library(reshape)
library(lattice)
library(grid)
library(gridBase)
source("functions/color-dendrogram.r")
source("functions/classify.r")

source("userinputs.r")

Arg <- tail(commandArgs(),1)
runno <- as.integer(Arg)
if( is.na(runno) )
  stop("enter run number as integer")

trymkdir <- function(x)
  if(!file.exists(x)) dir.create(x)

clusterout <- file.path(FOLDER,"Cluster")
clusterpath <- file.path(FOLDER,"Allplots","cluster")
trymkdir(dirname(clusterpath))
trymkdir(clusterpath)
trymkdir(clusterout)

dendpanel <- function(x, y, ..., dendrogram) {
  pushViewport(viewport(gp=gpar(fontsize=8)),
               viewport(y=unit(0.95, "npc"), width=0.9,
                        height=unit(0.8, "npc"),##unit(0.95, "npc") - space,
                        just="top"))
  par(plt=gridPLT(), new=TRUE, ps=8)
  plot(dendrogram)
  popViewport(2)
}

specpanel <- function(x,y,..., samples, groups,subscripts) {
  i <- as.character(groups[subscripts][1])
  col <- trellis.par.get("superpose.line")$col[i]
  for( elem in split(data.frame(x,y), samples[subscripts]) )
    with(elem,panel.lines(x,y,col=col,...))
}

normalize <- function(x)
  x/sqrt(crossprod(x,x))

###_2. load data

datamat <- read.delim(file.path(runnum(runno),"ExplainedVariation.txt"),
                      comment.char="",check.names=FALSE)
datamat$Resid <- NULL

###_3. clustering

###_ a. prepare matrix

##datamat[] <- t(apply(datamat,1,normalize))

###_ b. do the clustering

#ngrps <- max(nFactors)
hc <- hclust(dist(datamat))
## plot(hc)

for( ngrps in nFactors ) {
  print(ngrps)

  grps <- cutree(hc,k = ngrps)
  write.table(data.frame(sample=names(grps),cluster=grps),
              file.path(clusterout,sprintf("%s-%d-pmf-clusters.txt",basename(runnum(runno)),ngrps)),
              sep="\t",row.names=FALSE,quote=FALSE)

###_ c. plot it

###_  i. tree

  mycolors <- structure(rainbow(max(grps)),names=seq_len(max(grps)))
  p <- dendrapply(as.dendrogram(hc),colLab,grps,mycolors)
  ## plot(p)

  treeplot <- xyplot(y~x,data=data.frame(x=0,y=0),panel=dendpanel,
                     dendrogram=p,
                     par.settings=list(axis.line=list(col=NA)),
                     scales=list(draw=FALSE),xlab="",ylab="")

###_  ii. spectra
  
  long <- melt(as.matrix(datamat)*1e2,varnames=c("samples","variables"))
  long$grps <- factor(grps[long$samples])

  panel <- function(...,groups,subscripts,groupcols,numpts) {
    npts <- numpts[groups[subscripts]]
    color <- groupcols[groups[subscripts]]
    trellis.par.set(box.rectangle=list(col=color),
                    box.umbrella=list(col=color),
                    plot.symbol=list(col=color,cex=.5))
    panel.bwplot(...,col=color)
    grid.text(sprintf("n = %d",npts),x=1,y=1,hjust=1.1,vjust=1.5)
  }

  evplot <- bwplot(value ~ variables | grps, groups=grps, data=long,
                    samples=long$samples,jitter=TRUE,
                    panel=panel,
                    groupcols=mycolors,
                    numpts=table(grps),
                    as.table=TRUE,
                    scales=list(x=list(rot=30),
                      y=list(lim=c(0,100),at=seq(0,100,20))),
                    xlab="Factor",
                    ylab="Explained Variation (%)")
  
  ## combine
  grlayout <-
    grid.layout(1,2,widths=unit(c(.3,.7),"null"))

  trellis.device(png,file=file.path(clusterpath,
                       sprintf("%s-%d-pmf-clusters.png",basename(runnum(runno)),ngrps)),
                 width=10*96,height=5*96)
  grid.newpage()
  pushViewport(viewport(layout=grlayout))
  pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
  print(treeplot,newpage=FALSE)
  popViewport(1)
  pushViewport(viewport(layout.pos.col=2,layout.pos.row=1))
  print(evplot,newpage=FALSE)
  popViewport()
  dev.off()

}
