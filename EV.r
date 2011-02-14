## explained variance

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

###_* calculate explained variance
ExplainedVariance <- function(G,F,S,E) {
  num <- Reduce(`+`,Map(function(j)
                        abs(sweep(G,2,F[,j],`*`))/S[,j],
                        1:ncol(F)))
  res <- rowSums(abs(E)/S)
  den <- rowSums((Reduce(`+`,Map(function(h)
                                 abs(outer(G[,h],F[h,],`*`)),
                                 1:ncol(G)))+abs(E))/S)
  cbind(num,Resid=res)/den
}

###_* read, apply 
calcEV <- function(runno,FOLDER,export=TRUE) {
  ## read matrices
  dn <- list(samples=readLines(file.path(FOLDER,"samples.txt")),
             wavenumbers=readLines(file.path(FOLDER,"wavenumbers.txt")))
  mat <- as.matrix(read.table(file.path(FOLDER,"matrix.dat")))
  stdev <- as.matrix(read.table(file.path(FOLDER,"std_dev.dat")))
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
  evg <- ExplainedVariance(soln$G,soln$F,stdev,eij)

  if(export)
    write.table(evg,file=file.path(FOLDER,sprintf("%s_%s",basename(FOLDER),runno),
                  "ExplainedVariance.txt"),
            sep="\t",quote=FALSE)
  evg
}

runno <- sub(patt,"\\1",list.files(FOLDER,patt))
invisible(lapply(runno,calcEV,FOLDER=FOLDER))
