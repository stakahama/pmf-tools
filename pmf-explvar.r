####################
## PMF execution and postprocessing program
## ~EV.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


## explained variance

###_* command-line arguments
Arg <- tail(commandArgs(),1)
if( Arg=="--args" || Arg=="1" ) {
  FILENAME <- "ExplainedVariation.txt"
  FN <- quote(ExplainedVariation)
} else {
  FILENAME <- "ExplainedVariation_unnormalized.txt"
  FN <- quote(ExplainedVariation_unnormalized)
}

###_* import
source("userinputs.r")
source("filepaths.r")
source("functions/harvest.r")
patt <- ".+\\_([0-9]{3})"

###_* volumes
divbyv <- 
  if(is.null(projectinfofile)) identity else {
    local({
      pjinfo <- read.delim(projectinfofile)
      names(pjinfo) <- tolower(names(pjinfo))
      vols <- with(pjinfo,`names<-`(volume,filterid))
      function(m) m/vols[rownames(m)]
    })
  }

###_* calculate explained variation
ExplainedVariation <- function(G,F,S,E) {
  num <- Reduce(`+`,Map(function(j)
                        abs(sweep(G,2,F[,j],`*`))/S[,j],
                        1:ncol(F)))
  res <- rowSums(abs(E)/S)
  den <- rowSums((Reduce(`+`,Map(function(h)
                                 abs(outer(G[,h],F[h,],`*`)),
                                 1:ncol(G)))+abs(E))/S)
  cbind(num,Resid=res)/den
}

ExplainedVariation_unnormalized <- function(G,F,S,E) {
  num <- Reduce(`+`,Map(function(j)
                        abs(sweep(G,2,F[,j],`*`))/S[,j],
                        1:ncol(F)))
  res <- rowSums(abs(E))
  den <- rowSums((Reduce(`+`,Map(function(h)
                                 abs(outer(G[,h],F[h,],`*`)),
                                 1:ncol(G)))+abs(E)))
  cbind(num,Resid=res)/den
}

ExplainedVariation <- eval(FN)

###_* read, apply 
calcEV <- function(runno,FOLDER,export=TRUE) {
  ## read matrices
  dn <- list(samples=readLines(file.path(FOLDER,"samples.txt")),
             variables=readLines(file.path(FOLDER,"variables.txt")))
  mat <- as.matrix(read.table(file.path(FOLDER,"MATRIX.DAT")))
  stdev <- as.matrix(read.table(file.path(FOLDER,"STD_DEV.DAT")))
  dimnames(mat) <- dn
  dimnames(stdev) <- dn

  ## read solution
  soln <- getsoln(FOLDER, runno)

  ## divide by volume (this does not actually make a difference on the explained variance)
  mat[] <- divbyv(mat)
  soln$G[] <- divbyv(soln$G)

  ## compute residuals
  eij <- mat - with(soln,G %*% F)

  ## calculate explained variance
  evg <- ExplainedVariation(soln$G,soln$F,stdev,eij)

  if(export)
    write.table(evg,file=file.path(FOLDER,sprintf("%s_%s",basename(FOLDER),runno),
                  FILENAME),
            sep="\t",quote=FALSE)
  evg
}

runno <- sub(patt,"\\1",list.files(FOLDER,patt))

for( x in runno ) {
  if(!newsim &&
     file.exists(file.path(FOLDER,sprintf("%s_%s",basename(FOLDER),x),
                           FILENAME))) next()
  tryCatch({
    calcEV(x,FOLDER)
    print(x)
  },error=function(e) {
    print(sprintf("Error %s",x))
    print(e)
  })
}
