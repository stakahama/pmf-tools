populate.env("mylib","~/lib/R")
source("userinputs.r")
source("acsm-userinputs.r")

## load
here <- getwd()
setwd(file.path(Sys.getenv("HOME"),"VirtualBoxDocuments/ACSMdir"))  
load.itx(pmfmatfile)
org_specs <- t(org_specs)
orgspecs_err <- t(orgspecs_err)
acsm_time <- IGORdatetime(acsm_time)
ErrorMatB4DownWeight <- t(ErrorMatB4DownWeight)
setwd(here)

## subset functions
validsamples <- function() {
  ## local assignment
  between <- function(x,b) x >= min(b) & x < max(b)
  samples <- with(list(x=extendDiff(acsm_time)*24*60,
                       bounds=list(c(17,19),c(29,37))),
                  sapply(x,between,bounds[[1]]) | sapply(x,between,bounds[[2]]))
  samples <- samples & apply(org_specs[,amus > 150],1,function(x) all(x==0))
  ## global assignment
  org_specs <<- org_specs[samples,]
  orgspecs_err <<- orgspecs_err[samples,]
  acsm_time <<- acsm_time[samples]
  ErrorMatB4DownWeight <<- ErrorMatB4DownWeight[samples,]
}
validamus <- function() {
  ## local assignment
  minval <- 1e-5 ##.Machine$double.eps
  ## columns <- apply(ErrorMatB4DownWeight,2,function(x,m) any(x > m),minval)
  columns <- apply(orgspecs_err,2,function(x,m) any(x > m),minval)  
  ## global assignment  
  org_specs <<- org_specs[,columns]
  orgspecs_err <<- orgspecs_err[,columns]
  acsm_time <<- acsm_time[columns]
  ErrorMatB4DownWeight <<- ErrorMatB4DownWeight[,columns]
}
sdvgt0 <- function() {
  ## local assignment
  minval <- 1e-5 ##.Machine$double.eps
  ## samples <- apply(ErrorMatB4DownWeight,1,function(x,m) all(x > m),minval)
  samples <- apply(orgspecs_err,1,function(x,m) all(x > m),minval)
  ## global assignment    
  org_specs <<- org_specs[samples,]
  orgspecs_err <<- orgspecs_err[samples,]
  acsm_time <<- acsm_time[samples]
  ErrorMatB4DownWeight <<- ErrorMatB4DownWeight[samples,]
}
specialcase <- function() { ## --- june only ---
  samples <- acsm_time >= "06/01/10" ##& org_specs[,amus==44] > -.15
  select <- amus > min(amus)
  ## global assignment    
  org_specs <<- org_specs[samples,select]
  orgspecs_err <<- orgspecs_err[samples,select]
  acsm_time <<- acsm_time[samples]
  ErrorMatB4DownWeight <<- ErrorMatB4DownWeight[samples,select]
  amus <<- amus[select]
}
june_nonneg <- function() { ## --- june only ---
  samples <- acsm_time >= "06/01/10" ##& org_specs[,amus==44] > -.15
  select <- TRUE
  ## global assignment    
  org_specs <<- abs(org_specs[samples,select])
  orgspecs_err <<- orgspecs_err[samples,select]
  acsm_time <<- acsm_time[samples]
  ErrorMatB4DownWeight <<- ErrorMatB4DownWeight[samples,select]
  amus <<- amus[select]
}
june_lowmz <- function() { ## --- june only ---
  samples <- acsm_time >= "06/01/10" ##& org_specs[,amus==44] > -.15
  select <- amus < 100
  ## global assignment    
  org_specs <<- org_specs[samples,select]
  orgspecs_err <<- orgspecs_err[samples,select]
  acsm_time <<- acsm_time[samples]
  ErrorMatB4DownWeight <<- ErrorMatB4DownWeight[samples,select]
  amus <<- amus[select]
}
june_abslowmz <- function() { ## --- june only ---
  samples <- acsm_time >= "06/01/10" ##& org_specs[,amus==44] > -.15
  select <- amus < 100
  ## global assignment    
  org_specs <<- abs(org_specs[samples,select])
  orgspecs_err <<- orgspecs_err[samples,select]
  acsm_time <<- acsm_time[samples]
  ErrorMatB4DownWeight <<- ErrorMatB4DownWeight[samples,select]
  amus <<- amus[select]
}
mzle100 <- function() { ## --- june only ---
  select <- amus <= 100
  ## global assignment    
  org_specs <<- abs(org_specs[,select,drop=FALSE])
  orgspecs_err <<- orgspecs_err[,select,drop=FALSE]
  amus <<- amus[select]
  ErrorMatB4DownWeight <<- ErrorMatB4DownWeight[,select,drop=FALSE]
}
mzle125 <- function() { ## --- june only ---
  select <- amus <= 125
  ## global assignment    
  org_specs <<- abs(org_specs[,select,drop=FALSE])
  orgspecs_err <<- orgspecs_err[,select,drop=FALSE]
  amus <<- amus[select]
  ErrorMatB4DownWeight <<- ErrorMatB4DownWeight[,select,drop=FALSE]
}

## apply functions

## validsamples()
## sdvgt0()
## specialcase()
## june_nonneg()
## june_lowmz()
## june_abslowmz()
## mzle100()
mzle125()

## extract

if( !exists(errormatrixname) )
  errormatrixname <- "orgspecs_err"

runpath <- FOLDER
dir.create(runpath)
write.table(formatC(org_specs,format="g"),
            file=file.path(runpath,"MATRIX.DAT"),sep="\t",
            row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(formatC(eval(parse(text=errormatrixname)),format="g"),
            file=file.path(runpath,"STD_DEV.DAT"),sep="\t",
            row.names=FALSE,col.names=FALSE,quote=FALSE)
write(amus,file.path(runpath,"variables.txt"),ncol=1)
writeLines(format_chron(acsm_time,"%Y-%m-%d %T"),
           file.path(runpath,"samples.txt"))
