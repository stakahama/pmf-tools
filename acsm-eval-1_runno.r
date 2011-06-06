###_* library
library(lattice)
library(latticeExtra)
library(grid)
library(reshape)
library(chron)
library(fields)
source("functions/classify.r")
source("functions/stlib.r")

###_* data

source("userinputs.r")
source("acsm-userinputs.r")
OUTPATH <- file.path(FOLDER,"Allplots")

Arg <- tail(commandArgs(),1)
runno <- as.integer(Arg)
if( is.na(runno) )
  stop("enter run number as integer")

runno <- sprintf("%s_%03d",basename(FOLDER),runno)

###_ . sticks

print(load(orgsticksfile))
sticks2fmz <- function(sticks,RIE.Org=1.4,IE.NO3=2.3e-11)
  within(sticks,
         orgstickmatrix[] <- orgstickmatrix/RIE.Org/IE.NO3)
stickscorrection <- function(sticks)
  within(sticks,
         orgstickmatrix[] <- sweep(orgstickmatrix,1,usercorrection,`/`))
sticksaddrownames <- function(sticks)
  within(sticks,rownames(orgstickmatrix) <-
         format_chron(IGORdatetime(utctime)))
sticks2org <- function(amus,orgstickmatrix,select=10:125)
  rowSums(orgstickmatrix[,na.omit(match(select,amus))],na.rm=TRUE)
mzselector <- function(sticks) function(mzs)
  with(sticks,orgstickmatrix[,na.omit(match(mzs,amus)),drop=FALSE])
sticks[] <- sticksaddrownames(sticks)
orgsticks <- sticks
orgsticks[] <- stickscorrection(sticks2fmz(orgsticks))
selectmzs <- mzselector(orgsticks)
normalizeorg <- function(x)
  sweep(x,1,rowSums(selectmzs(eval(formals(sticks2org)$select))),`/`)


###_ . ACSM PMF

getpmfclass <- function(x) {
  m <- pmfmatches[x,"Class"]
  i <- which(m==m[duplicated(m)])
  if(length(i)==0) return(m)
  replace(m,i,sprintf("%s_%d",m[i],1:length(i)))
}

pmfmatches <- read.delim(file.path(FOLDER,"Allplots/PMF-matches.txt"),
                         row.names=1)
##
EV <- read.delim(file.path(FOLDER,runno,"ExplainedVariation.txt"),
                 check.names=FALSE)
EV[] <- sweep(EV,1,100,`*`)
##
pmfamus <- scan(file.path(FOLDER,"variables.txt"),0,quiet=TRUE)
G <- read.table(file.path(FOLDER,runno,"G_FACTOR.TXT"))
rownames(G) <- rownames(EV)
colnames(G) <- head(colnames(EV),-1)
F <- read.table(file.path(FOLDER,runno,"F_FACTOR.TXT"))
colnames(F) <- sprintf("mz%d",pmfamus)
rownames(F) <- colnames(G)
pmforg <- sapply(rownames(F),function(x,mz,G,F)
                 sticks2org(mz,matrix(G[,x]) %*% t(F[x,])),
                 pmfamus,as.matrix(G),as.matrix(F))
rownames(pmforg) <- rownames(G)
##
###_* merging function
combine <- function(TS,MZ) {
  TS <- as.data.frame(TS)
  MX <- as.data.frame(MZ)
  key <- rownames(TS)
  merge(rename(melt(cbind(key=key,TS),
                    id=1),c(variable="factor",value="strength")),
        rename(melt(cbind(key=key,as.data.frame(MZ[key,])),
                    id=1),c(variable="mz",value="signal")),
        "key")
}

###_* plotting function
usepts <- function(x)
  findInterval(x,quantile(x,c(.02,.98),na.rm=TRUE))==1

prepanel <- function(x,y,...) {
  i <- usepts(x)
  x <- x[i]; y <- y[i]
  list(xlim=range(x,na.rm=TRUE),ylim=range(y,na.rm=TRUE))
}

panel <- function(x,y,...) {
  i <- usepts(x)  
  x <- x[i]; y <- y[i]
  panel.xyplot(x,y,cex=.5,...)
  ## panel.smoothScatter(x,y,...)
  grid.text(sprintf("r = %.2f",cor(x,y,use="pairwise")),
            1,0,just=c(1,0),gp=gpar(fontface="italic"))
}

hplot <- function(bf,...,col=c("darkblue","darkred"),e=.15) {
  x <- as.numeric(substring(colnames(bf),3))
  plot(outer(x,0.5*nrow(bf)*seq(-1,1,,nrow(bf))*e,`+`),t(bf),type="h",lwd=3,lend=3,
       col=rep(col,each=length(x)),
       xlab="m/z",ylab="Correlation (r)",...)
}

tspanel1 <- function(x,y,...,groups,subscripts) {
  qy <- sapply(split(y,x),quantile,c(.25,.5,.75),na.rm=TRUE)
  newx <- as.numeric(levels(x))
  panel.points(newx,qy[2,])
  panel.segments(newx,qy[1,],newx,qy[3,])
}

tsprepanel2 <- function(.p) function(x,y,...) {
  list(ylim=c(0,quantile(y,.p,na.rm=TRUE)))
}

tspanel2 <- function(x,y,...,groups,subscripts) {
  for( iter in zip(head(trellis.par.get("superpose.symbol")$col,nlevels(groups)),
                   split(data.frame(x,y),groups[subscripts])) )
    with(iter[[2]],{
      qy <- sapply(split(y,x),quantile,c(.25,.5,.75),na.rm=TRUE)
      newx <- as.numeric(levels(x))
      panel.points(newx,qy[2,],col=iter[[1]])
      panel.segments(newx,qy[1,],newx,qy[3,],col=iter[[1]])
    })
}

###_* plots

###_ . scatter

mg <- combine(pmforg,selectmzs(c(43,44,57,60)))
levels(mg$factor) <- getpmfclass(levels(mg$factor))
sout <- xyplot(strength ~ signal | mz*factor, data=mg,
              prepanel=prepanel,panel=panel,
              xlab=expression("m/z Org"~(mu*g/m^3)),
              ylab=expression("PMF Org"~(mu*g/m^3)),
              scales=list(relation="free"))

pdf(file.path(OUTPATH,sprintf("%s-scatter.pdf",runno)),width=8,height=6)
print(useOuterStrips(combineLimits(sout)))
dev.off()

###_ . correlations

pdf(file.path(OUTPATH,sprintf("%s-correlations.pdf",runno)),width=8,height=5)
par(mfrow=c(1,1),mar=c(4,4,1,1),mgp=c(2.2,.5,0))
bf <- cor(pmforg,selectmzs(12:100)[rownames(pmforg),],use="pairwise")
hplot(bf,ylim=c(-1,1),xaxt="n",col=tim.colors(nrow(bf)))
axis(1,seq(12,100,4))
abline(h=0,col=8)
legend("bottomright",fill=tim.colors(nrow(bf)),getpmfclass(colnames(pmforg)))
dev.off()

###_ . timeseries

###_  : absolute
orgts <- with(list(psttime=as.chron(rownames(pmforg))-8/24),
              append(data.frame(pmforg,check.names=FALSE),
                     data.frame(psttime=psttime,
                                hour=factor(hours(psttime)),
                                weekdays=weekdays(psttime)),0))
lg <- melt(orgts,id=1:3)
levels(lg$variable) <- getpmfclass(levels(lg$variable))
xout <- xyplot(value ~ hour, data=lg,
               par.settings=list(superpose.symbol=
                 list(col=tim.colors(nlevels(lg$variable)))),
               groups = variable,
               prepanel=tsprepanel2(.98),
               panel=tspanel2,
               auto=list(space="right"),
               xlab="Hour (PST)",ylab=expression("Org"~(mu*g/m^3)))

###_  : fraction
orgfrac <- funcall(function(x) data.frame(hour=factor(as.numeric(rownames(x))),
                                          sweep(x,1,rowSums(x),`/`),
                                          check.names=FALSE),
                   t(sapply(split(as.data.frame(pmforg),orgts$hour),colSums)))
lgn <- melt(orgfrac,id=1)
levels(lgn$variable) <- getpmfclass(levels(lgn$variable))
bout <- barchart(value~hour,data=lgn,groups=variable,
         par.settings=list(superpose.polygon=
           list(col=tim.colors(nlevels(lgn$variable)))),
         stack=TRUE,auto=list(space="right",reverse.rows=TRUE),
         scales=list(y=list(lim=c(0,1),at=seq(0,1,.2))),
         box.width=1,
         xlab="Hour (PST)",ylab="Org fraction")

pdf(file.path(OUTPATH,sprintf("%s-tseries.pdf",runno)),width=8,height=5)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1)))
pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
print(xout,newpage=FALSE)
popViewport()
pushViewport(viewport(layout.pos.col=1,layout.pos.row=2))
print(bout,newpage=FALSE)
popViewport()
dev.off()
