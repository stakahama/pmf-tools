###_* library
library(RSQLite)
source("functions/classify.r")
source("userinputs.r")
## inputs: FOLDER, dbname

if(!exists("dbname")) dbname <- "dbfiles/acsm-refspec.db"

###_* load reference spectra
db2matrix <- function(refspec)
  `rownames<-`(as.matrix(refspec[,-(1:2)]),
               do.call(paste,c(refspec[,(1:2)],list(sep=":"))))


tablename <- "reference_spectra"
drv <- dbDriver("SQLite")
conn <- dbConnect(drv,dbname=dbname)
refspec <- dbReadTable(conn,tablename)
dbDisconnect(conn)

refspec <- db2matrix(refspec)
dowant <- c("OOA","HOA","BB")
refspec <-
  refspec[rownames(refspec) %in%
          unlist(lapply(dowant,grep,
                        rownames(refspec),
                        value=TRUE,ignore.case=TRUE)),]
refvars <- as.integer(sub("^mz","",colnames(refspec)))

###_* inputs

## labels
label <- c(18,27,28,29,32,43,44,55,57,60,73,80,94)
refxval <- list(x=refvars[refvars %in% label],
                cols=refvars %in% label)
## select m/zs
## s <- refvars < 100
s <- TRUE ## (for plotting)
select <- sprintf("mz%d",c(43,44,57,60,73))

simgrid <- read.delim(file.path(FOLDER,"simgrid.txt"),row.names=1)
simgrid <- simgrid[with(simgrid,order(nFactors,FPEAK,Seed)),]
simgrid <- unique(simgrid)
sims <- rownames(subset(simgrid,nFactors %in% 2:8))

correlations <- structure(vector("list",length(sims)),names=sims)
dir.create(file.path(FOLDER,"Allplots"))
pdf(file.path(FOLDER,"Allplots",sprintf("PMF-matches_%s.pdf",basename(FOLDER))))
for( g in sims ) {
  runno <- as.integer(sub(".+\\_([0-9])","\\1",g))
  X <- as.matrix(read.table(file.path(runnum(runno),"F_FACTOR.TXT")))
  samples <- sprintf("%s-%02d",basename(runnum(runno)),1:nrow(X))
  vars <- scan(file.path(FOLDER,"variables.txt"),0,quiet=TRUE)
  dimnames(X) <- list(samples,sprintf("mz%d",vars))
  xval <- list(x=vars[vars %in% label],
               cols=vars %in% label)
  ##
  ## select <- colnames(X)
  ## select is now defined outside of loop (fixed set of values)

  ref.X <-
    outer(structure(1:nrow(refspec),names=rownames(refspec)),
          structure(1:nrow(X),names=rownames(X)),
          function(i,j,X,Y)
          mapply(function(.i,.j) corr(p2_norm(X[.i,]),p2_norm(Y[.j,])),i,j),
          refspec[,select],X[,select])
  ##
  near <- structure(rownames(ref.X)[apply(ref.X,2,which.max)],
                    names=colnames(ref.X))
  correlations[[g]] <- Map(function(a,b) structure(ref.X[a,b],names=a),near,names(near))

  par(mfrow=c(length(near),2),mar=c(0,0,0,0),oma=c(4,4,2,2),mgp=c(2.2,.5,0))
  for( i in 1:length(near) ) {
    plot(vars,X[names(near)[i],],type="h",ylim=c(0,1.2*max(X[names(near)[i],])),
         xaxt="n",xlim=c(8,100),yaxs="i",col=2)
    text(xval[[1]],X[i,xval[[2]]],adj=c(.5,0),xval[[1]])
    text(par("usr")[2],par("usr")[4],adj=c(1,1),rownames(X)[i],cex=1.2)  
    axis(1,,if( i==length(near) ) TRUE else FALSE)
    plot(refvars[s],refspec[near[i],s],type="h",ylim=c(0,1.2*max(refspec[near[i],],na.rm=TRUE)),
         xaxt="n",xlim=c(8,100),yaxt="n",yaxs="i",col=4)
    axis(4)
    text(refxval[[1]],refspec[near[i],refxval[[2]]],adj=c(.5,0),refxval[[1]])
    text(par("usr")[2],par("usr")[4],adj=c(1,1),near[i],cex=1.2)
    text(par("usr")[2],par("usr")[4]-par("cxy")[2],adj=c(1,1),cex=1.2,
         sprintf("r=%.2f",ref.X[near[i],names(near)[i]]))
    axis(1,,if( i==length(near) ) TRUE else FALSE)  
  }
  mtext("m/z",1,outer=TRUE,line=2)
  mtext("Signal",2,outer=TRUE,line=2)
  mtext(do.call(paste,c(as.list(sprintf("%s=%s",names(simgrid),simgrid[g,])),sep=",")),
        3,outer=TRUE)
}
dev.off()

###_* summary plot
library(lattice)
library(latticeExtra)
library(reshape)
library(RColorBrewer)

OAclass <- function(x,ref=rownames(get("refspec",environment(OAclass)))) {
  fn <- function(x) 
    sapply(strsplit(x,"[:]"),function(.x)
           sub(".*((H|O|BB)OA(1|2|I|II)?).*","\\1",.x[2]))
  sb <- function(x) Reduce(function(x,y) sub(y[1],y[2],x),
                           list(c("II","2"),c("I","1")),x)
  factor(sb(fn(x)),unique(sb(fn(ref))))
}

dfr <- do.call(rbind,Map(function(x,y) 
                         data.frame(sim=y,factor=names(x),
                                    r=c(unlist(x)),
                                    match=sapply(x,names)),
                         correlations,names(correlations)))
rownames(dfr) <- NULL
dfr <- cbind(dfr,simgrid[dfr$sim,])
dfr$Class <- OAclass(dfr$match)
dfr <- within(dfr,{
  nFactors <- factor(nFactors)
  FPEAK <- factor(FPEAK)  
  Seed <- factor(Seed)
})

##
out <- xyplot(r~FPEAK | nFactors*Seed, groups=Class,data=dfr,
              ylim=c(0.45,1.05),
              par.settings=list(superpose.symbol=
                list(pch=1:length(levels(dfr$Class)))),
              scales=list(x=list(rot=90),
                y=list(at=seq(0,1,.1))),
              auto.key=list(space="right"),#list(x=.7,y=.3),
              as.table=TRUE)
pdf(file.path(FOLDER,"Allplots",sprintf("summaryplot_%s.pdf",basename(FOLDER))),
    width=10,height=6)
useOuterStrips(out)
dev.off()

##
wide <- recast(dfr[,c("sim",names(simgrid),"Class")],
               as.formula(paste(paste(c("sim",names(simgrid)),collapse="+"),
                                "~ value")),
               measure.var="Class")
class(wide) <- "data.frame"
##
long <- melt(wide,id.var=c("sim",names(simgrid)),variable_name="Class")

mycol <- brewer.pal(max(long$value)+1,"Blues")
## mycol[1] <- rgb(1,1,1)
out <- levelplot(value ~ FPEAK+Class | nFactors*Seed, data=long,
                 scales=list(x=list(rot=90)),
                 at=-1:max(long$value)+.5,
                 col.regions=mycol,
                 colorkey=list(labels=list(at=0:max(long$value))),
                 as.table=TRUE)
pdf(file.path(FOLDER,"Allplots",sprintf("heatmap_%s.pdf",basename(FOLDER))),
    width=10,height=6)
useOuterStrips(out)
dev.off()
