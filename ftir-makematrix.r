####################
## PMF execution and postprocessing program
## ~ftir-makematrix.r~
## $Rev: 7 $
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################

source("ftir-inputs.r")
dir.create(datapath)

source("functions/cleaning.r")

load(file.path(ftirpath,outputpath,sprintf("%s-bl.rda",projectname)))
blanks <- read.csv(file.path(ftirpath,outputpath,"blanksummary.csv"),
                   row.names=1)

for( i in 1:nrow(bl) ) 
  bl[i,c("wavenumber","baselined")] <-
  with(bl[i,],cleanbl(wavenumber,baselined))

##{{{ --- projectinfo file ---
## pjinfo <- read.delim(file.path(ftirpath,inputpath,pjfile),as.is=TRUE)
## common <- intersect(rownames(bl),pjinfo$FilterID)
## pjinfo <- pjinfo[match(common,pjinfo$FilterID),]
## bl <- bl[common,]
## samples.use <- with(pjinfo,FilterID[Type=="S"])
##}}}

## Subset Type == "S"
## Make sample matrix
## samples.use should be defined in ftirinputs.r
samples.use <- if(is.character(samples.use) && file.exists(samples.use))
  readLines(samples.use) else samples.use
mat <- do.call(rbind,bl[samples.use,"baselined"])

## Select wavenumbers
deps <- with(.Machine,max(abs(c(double.eps,double.neg.eps))))
wavenums.use <- with(blanks,Sd > deps)
mat <- mat[,wavenums.use]
stdev <- outer(rep(1,nrow(mat)),blanks$Sd[wavenums.use])
## stdev <- 1/mat

form <- function(x) structure(sprintf("%.5f",x),dim=dim(x))

write.table(form(mat),file.path(datapath,"matrix.dat"),
            row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
write.table(form(stdev),file.path(datapath,"std_dev.dat"),
            row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
writeLines(format(bl[[1,"wavenumber"]][wavenums.use],nsmall=3),file.path(datapath,"variables.txt"))
writeLines(rownames(mat),file.path(datapath,"samples.txt"))
