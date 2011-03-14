####################
## PMF execution and postprocessing program
## ~harvest.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


## harvest PMF matrices

funcall <- function(FUN,...) FUN(...)

getsoln <- function(FOLDER,solution,xvariables) {
  path <- unique(file.path(FOLDER,sprintf("%s_%s",basename(FOLDER),solution)))
  g <- funcall(function(p) {
    samples <- readLines(file.path(dirname(p),"samples.txt"))
    g <- as.matrix(read.table(file.path(p,"G_FACTOR.TXT")))
    `dimnames<-`(g,list(samples,
                        sprintf("%s-%02d",basename(p),1:ncol(g))))
  },path)
  f <- funcall(function(p) {
    wn <- readLines(file.path(dirname(p),xvariables))
    f <- matrix(scan(file.path(p,"F_FACTOR.TXT"),quiet=TRUE),ncol=length(wn),
                byrow=TRUE)
    `dimnames<-`(f,list(sprintf("%s-%02d",basename(p),1:nrow(f)),wn))
  },path)
  list(F=f,G=g)
}
