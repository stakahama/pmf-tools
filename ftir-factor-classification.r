###_* library
source("functions/classify.r")

###_* load reference spectra
## --- please request ftir-refspec.db from s.t. and place it in a
## directory named dbfiles/ ---
library(RSQLite)
library(reshape)
loadspec <- function(specmat,vars) {
  wide <- cast(specmat,Spectrum~VariableName,value="Absorbance")
  structure(as.matrix(wide[,-1]),
            dimnames=list(wide[,1],
              vars$Wavenumber[match(names(wide)[-1],vars$VarName)]))
}
drv <- dbDriver("SQLite")
conn <- dbConnect(drv,dbname="dbfiles/ftir-refspec.db")
refspec <- loadspec(dbReadTable(conn,"spectramatrix"),
                    dbReadTable(conn,"wavenumbertable"))
refvars <- as.numeric(colnames(refspec))
refclasses <- dbReadTable(conn,"classes")
dbDisconnect(conn)

###_* inputs

attach(NULL,name="userinputs")
sys.source("~/VirtualBoxDocuments/programs/pmf/userinputs.r",
           as.environment("userinputs"))

## this is to select a subset of reference variables if desired
refxval <- list(x=refvars,cols=TRUE)
s <- TRUE

simgrid <- read.delim(file.path(FOLDER,"simgrid.txt"),row.names=1)
simgrid <- simgrid[with(simgrid,order(nFactors,FPEAK,Seed)),]
sims <- rownames(subset(simgrid,nFactors %in% 2:6))

align <- with(list(refspec=refspec,refvars=refvars), function(vars) {
  structure(t(apply(refspec,1,function(.x,r,v) approx(r,.x,v)$y,refvars,vars)),
            dimnames=list(rownames(refspec),sprintf("%.2f",vars)))
})

extend <- function(x,e=1.04) {
  .x <- na.omit(x)
  .e <- diff(range(.x))*(e-1)
  c(min(.x)-.e,max(.x)+.e)
}

dir.create(file.path(FOLDER,"Allplots")
correlations <- structure(vector("list",length(sims)),names=sims)
pdf(file.path(FOLDER,"Allplots","db-factor-matches.pdf"))
for( g in sims ) {
  runno <- as.integer(sub(".+\\_([0-9])","\\1",g))
  X <- as.matrix(read.table(file.path(runnum(runno),"F_FACTOR.TXT")))
  samples <- sprintf("%s-%02d",basename(runnum(runno)),1:nrow(X))
  vars <- scan(file.path(FOLDER,"variables.txt"),0,quiet=TRUE)
  ## dimnames(X) <- list(samples,sprintf("mz%d",vars))
  dimnames(X) <- list(samples,sprintf("%.2f",vars))
  ## xval <- list(x=vars[vars %in% label],
  ##              cols=vars %in% label)
  xval <- list(x=vars,cols=TRUE)
  ##
  rv <- align(vars)
  select <- (apply(rv,2,function(.x) all(!is.na(.x))) &
             apply(X,2,function(.x) all(!is.na(.x))))
  ## select <- sprintf("mz%d",c(29,43,44,57,60,73))
  ref.X <-
    outer(structure(1:nrow(refspec),names=rownames(refspec)),
          structure(1:nrow(X),names=rownames(X)),
          function(i,j,X,Y)
          mapply(function(.i,.j) corr(p2_norm(X[.i,]),p2_norm(Y[.j,])),i,j),
          rv[,select],X[,select])
  ## ---
  near <- structure(rownames(ref.X)[apply(ref.X,2,which.max)],
                    names=colnames(ref.X))
  correlations[[g]] <- Map(function(a,b) structure(ref.X[a,b],names=a),near,names(near))
  next()
  ##
  xlim <- c(4000,1400)
  par(mfrow=c(length(near),2),mar=c(0,0,0,0),oma=c(4,4,2,2),mgp=c(2.2,.5,0))
  for( i in 1:length(near) ) {
    plot(vars,X[names(near)[i],],type="l",ylim=extend(X[names(near)[i],]),
         xaxt="n",xlim=xlim,yaxs="i",col=2)
    text(par("usr")[2],par("usr")[4],adj=c(1,1),rownames(X)[i],cex=1.2)  
    axis(1,,if( i==length(near) ) TRUE else FALSE)
    plot(refvars[s],refspec[near[i],s],type="l",ylim=extend(refspec[near[i],]),
         xaxt="n",xlim=xlim,yaxt="n",yaxs="i",col=4)
    axis(4)
    text(par("usr")[2],par("usr")[4],adj=c(1,1),near[i],cex=1.2)
    text(par("usr")[2],par("usr")[4]-par("cxy")[2],adj=c(1,1),cex=1.2,
         sprintf("r=%.2f",ref.X[near[i],names(near)[i]]))
    axis(1,,if( i==length(near) ) TRUE else FALSE)  
  }
  mtext("Wavenumber",1,outer=TRUE,line=2)
  mtext("Absorbance",2,outer=TRUE,line=2)
  mtext(do.call(paste,c(as.list(sprintf("%s=%s",names(simgrid),simgrid[g,])),sep=",")),
        3,outer=TRUE)
}
dev.off()

###_* summary plot
library(lattice)
library(latticeExtra)
library(reshape)
library(RColorBrewer)

dfr <- do.call(rbind,Map(function(x,y) 
                         data.frame(sim=y,factor=names(x),
                                    r=c(unlist(x)),
                                    match=sapply(x,names)),
                         correlations,names(correlations)))
rownames(dfr) <- NULL
dfr <- cbind(dfr,simgrid[dfr$sim,])
dfr$Class <- factor(toupper(sub("(.+)\\-(.+)","\\1",dfr$match)),
                    refclasses$ID)
dfr <- within(dfr,{
  nFactors <- factor(nFactors)
  FPEAK <- factor(FPEAK)  
  Seed <- factor(Seed)
})

out <- xyplot(r~FPEAK | nFactors*Seed, groups=Class,data=dfr,
              ylim=c(0.45,1.05),
              par.settings=list(superpose.symbol=
                list(pch=1:length(levels(dfr$Class)))),
              scales=list(x=list(rot=90),
                y=list(at=seq(0,1,.1))),
              auto.key=list(space="right"),#list(x=.7,y=.3),
              as.table=TRUE)
pdf(file.path(FOLDER,"Allplots","db-summaryplot.pdf"),width=8,height=8)
useOuterStrips(out)
dev.off()

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
pdf(file.path(FOLDER,"Allplots","db-heatmap.pdf"),width=8,height=8)
useOuterStrips(out)
dev.off()
