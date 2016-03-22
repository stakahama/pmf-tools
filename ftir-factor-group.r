###_* libraries and functions

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

## library(reshape)
library(latticeExtra)
library(lattice)

specpanel <- function(x,y,..., samples, groups,subscripts) {
  i <- as.character(groups[subscripts][1])
  col <- trellis.par.get("superpose.line")$col[i]
  for( elem in split(data.frame(x,y), samples[subscripts]) )
    with(elem,panel.lines(x,y,col=col,...))
}

read.F <- function(x) {
  run <- basename(dirname(x))
  ff <- as.matrix(read.table(x))
  rownames(ff) <- sprintf("%s-%02d",run,1:nrow(ff))
  ff
}

normalize <- function(x)
  x/sqrt(crossprod(x,x))

###_* inputs

fi <- list.files(FOLDER,"F_FACTOR.TXT",rec=TRUE,full=TRUE)
mat <- do.call(rbind,Map(read.F,fi))
colnames(mat) <- readLines(file.path(FOLDER,"variables.txt"))

###_* normalize, cluster

mat[] <- t(apply(mat,1,normalize))
hc <- hclust(dist(mat))

###_* trim tree, long data form

ngrps <- max(nFactors)
grp <- cutree(hc,k=ngrps)
numFac <- with(list(x=sub("-[0-9]+$","",rownames(mat))), table(x)[x])
longspec <- data.frame(numFac=factor(rep(numFac,ncol(mat))),
                       grp=factor(rep(grp,ncol(mat))),
                       comp=rep(rownames(mat),ncol(mat)),
                       variable=rep(as.numeric(colnames(mat)),each=nrow(mat)),
                       value=c(mat))

###_* plot

mycolors <- with(list(n=nlevels(longspec$grp)),
                 setNames(rainbow(n),seq(n)))
specplot <- xyplot(value ~ variable | grp*numFac, groups=grp, data=longspec,
                   panel=specpanel, samples=longspec$comp,
                   type="l", as.table=TRUE,
                   par.settings=list(superpose.line=list(col=mycolors)),
                   xlim=c(4000,1000),
                   xlab=expression("Wavenumber"~(cm^-1)),
                   ylab="Absorbance")

pdf(file.path(FOLDER,"Allplots",sprintf("factors-%d-clusters.pdf",ngrps)),
    width=8,height=5)
print(useOuterStrips(specplot))
dev.off()

###_* export

write.table(convertFactor(longspec[,c("comp","grp")], 2),
            file.path(FOLDER, "Allplots", sprintf("factors-%d-clusters.txt", ngrps)),
            sep="\t", row.names=FALSE)


