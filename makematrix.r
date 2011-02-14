## setwd("~/Desktop/pmf_nt")
ftirpath <- "../CalNex/ftir"
pjname <- "CalNex"
outputpath <- "outputs"
datapath <- "runs/CalNex1"
dir.create(datapath)

source("functions/cleaning.r")
load(file.path(ftirpath,outputpath,sprintf("%s-bl.rda",pjname)))
blanks <- read.csv(file.path(ftirpath,outputpath,"blanksummary.csv"),row.names=1)

for( i in 1:nrow(bl) ) 
  bl[i,c("wavenumber","baselined")] <-
  with(bl[i,],cleanbl(wavenumber,baselined))

## Subset Type == "S"
## Make sample matrix
samples.use <- TRUE
mat <- do.call(rbind,bl[samples.use,"baselined"])

## Select wavenumbers
deps <- with(.Machine,max(abs(c(double.eps,double.neg.eps))))
wavenums.use <- with(blanks,Sd > deps)
## wavenums.use <- apply(mat > deps),2,all)
mat <- mat[,wavenums.use]
stdev <- outer(rep(1,nrow(mat)),blanks$Sd[wavenums.use])
## stdev <- 1/mat

form <- function(x) structure(sprintf("%.5f",x),dim=dim(x))

write.table(form(mat),file.path(datapath,"matrix.dat"),
            row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
write.table(form(stdev),file.path(datapath,"std_dev.dat"),
            row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
writeLines(format(bl[[1,"wavenumber"]][wavenums.use],nsmall=3),file.path(datapath,"wavenumbers.txt"))
writeLines(rownames(mat),file.path(datapath,"samples.txt"))
