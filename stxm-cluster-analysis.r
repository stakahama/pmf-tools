###_1. load libraries
library(reshape)
library(lattice)
library(grid)
library(gridBase)
source("functions/color-dendrogram.r")

source("userinputs.r")

trymkdir <- function(x)
  if(!file.exists(x)) dir.create(x)

clusterpath <- file.path(FOLDER,"Allplots","cluster")
trymkdir(dirname(clusterpath))
trymkdir(clusterpath)

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

edgenorm <- function(x,energies) {
  preedge <- mean(x[energies < 283],na.rm=TRUE)
  postedge <- mean(x[energies > 305]-preedge,na.rm=TRUE)
  normalize((x-preedge)/postedge)
}

###_2. load data
datamat <-
  as.matrix(read.table(file.path(FOLDER,"matrix.dat"),colClasses="numeric"))
dimnames(datamat) <-
  list(samples=readLines(file.path(FOLDER,"samples.txt")),
       variables=readLines(file.path(FOLDER,"variables.txt")))

###_3. clustering

###_ a. prepare matrix

datamat[] <- t(apply(datamat,1,edgenorm,
                     as.numeric(dimnames(datamat)$variables)))

###_ b. do the clustering

#ngrps <- max(nFactors)
hc <- hclust(dist(datamat))
## plot(hc)

for( ngrps in nFactors ) {
  print(ngrps)
  
  grps <- cutree(hc,k = ngrps)

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
  
  longspec <- melt(datamat)
  longspec$grps <- factor(grps[longspec$samples])

  specplot <- xyplot(value ~ variables | grps, groups=grps, data=longspec,
                     panel=specpanel, samples=longspec$samples,
                     type="l", as.table=TRUE,
                     par.settings=list(superpose.line=list(col=mycolors)),
                     xlab="Energy (eV)",
                     ylab="Absorbance")
  
###_  iii. image
  
  pixels <- data.frame(structure(t(mapply(as.numeric,strsplit(dimnames(datamat)$samples,","))),dimnames=list(NULL,c("x","y"))),grps=grps)

  implot <- levelplot(grps~x+y,data=pixels,
                      at=0:length(mycolors)+.5,
                      colorkey=list(labels=list(at=1:length(mycolors),
                                      lab=1:length(mycolors))),
                      par.settings=list(regions=list(col=mycolors)),
                      xlab=expression("x"~(mu*m)),
                      ylab=expression("x"~(mu*m)))

  ## combine
  grlayout <-
    grid.layout(2,2,heights=unit(rep(.5, 2),"null"),
                widths=unit(rep(.5, 2),"null"))

  trellis.device(png,file=file.path(clusterpath,
                       sprintf("%d-clusters.png",ngrps)),
                 width=8*96,height=8*96)
  grid.newpage()
  pushViewport(viewport(layout=grlayout))
  pushViewport(viewport(layout.pos.col=1,layout.pos.row=2))
  print(treeplot,newpage=FALSE)
  popViewport(1)
  pushViewport(viewport(layout.pos.col=1:2,layout.pos.row=1))
  print(specplot,newpage=FALSE)
  popViewport()
  pushViewport(viewport(layout.pos.col=2,layout.pos.row=2))
  print(implot,newpage=FALSE)
  popViewport()
  dev.off()

}
