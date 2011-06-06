####################
## PMF execution and postprocessing program
## ~Rfuncs.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


## make simgrid
creategrid <- function(nFactors,FPEAK,Seeds,newsim) {
  gr <- unique(expand.grid(nFactors=nFactors,FPEAK=FPEAK,Seed=Seeds))
  if(newsim) {
    mx <- 0
  } else {
    if(file.exists(file.path(FOLDER,"simgrid.txt"))) {
      patt <- sprintf("%s\\_([0-9]+)",basename(FOLDER))
      simgrid <- read.delim(file.path(FOLDER,"simgrid.txt"),row.names=1,
                            colClasses=c("character","integer","numeric","numeric"))
      gr <- gr[!do.call(paste,gr) %in% do.call(paste,simgrid),]
      mx <- max(as.integer(sub(patt,"\\1",rownames(simgrid)))    )    
    } else {
      mx <- 0
    }
  }
  gr <- gr[with(gr,order(nFactors,FPEAK,Seed)),]  
  rownames(gr) <- paste(basename(FOLDER),
                        substring(1000+mx+(1:nrow(gr)),2),sep="_")
  gr
}

## modify INI file to adjust parameters
## nrow, ncol, nFactors, FPEAK

editrecs <- function(txt,fc,ncol,charperfield=20) { #charperfield <- 14
  j <- grep(txt,fc)
  if( (val <- ncol * charperfield) > 2000) 
    list(j,sub(" 2000",formatC(val,width=7,format="d"),fc[j]))
  else NULL
}
editoutformat <- function(fx,ncol) {
  ix <- grep("57[ ]+\".+\"",fx)
  list(ix,sub(",150\\(",sprintf(",%d(",ncol),fx[ix]))
}
editdimensions <- function(x,ncol,nrow,nfact) {
  j <- grep("Dimensions",x)+1  
  val <- Reduce(function(init,x) sub(x[[1]],x[[2]],init),
                list(list("4 ",paste(nfact," ",sep="")),
                     list("20",ncol),
                     list("40",nrow)),x[j])
  list(j,val)
}
repl1 <- function(x,nrow,ncol,nfact=4) {
  dimensions <- editdimensions(x,ncol,nrow,nfact)
  fieldlen <- Filter(Negate(is.null),
                     lapply(c('"MATRIX\\.DAT','"STD_DEV\\.DAT',
                              '"PMF3[2-4]\\.DAT','"MISC\\.DAT',
                              '"G\\_FACTOR\\.TXT',
                              '"F\\_FACTOR.TXT',
                              '"TEMP.TXT',
                              '"[$]\\.DAT'),
                            editrecs,x,ncol))
  outformat <- editoutformat(x,ncol)
  Reduce(function(init,x)
         if(is.null(x)) init else replace(init,x[[1]],x[[2]]),
         c(list(dimensions),fieldlen,list(outformat)),x)
}
repl2 <- function(x,fpeak=0) {
  j <- grep("FPEAK",x)+1
  val <- sub("0.00000",format(fpeak,nsmall=5),x[j])
  replace(x,j,val)
}
repl3 <- function(x,seedval=1) {
  j <- grep("Seed",x)+1
  val <- sub(" 1 ",sprintf(" %d ",seedval),x[j])
  replace(x,j,val)
}
makeFunc <- function(ini,mat) {
  rows <- dim(mat)[1]
  cols <- dim(mat)[2]
  function(nfact=4,fpeak=0,seed=1)
    repl3(repl2(repl1(ini,rows,cols,nfact),fpeak),seed)
}
## goodstart, sortfactorsg, and "lims"
## sortfactorsf not used?
repOpts <- function(x,i) {
  optstring <- ifelse(i > 1,"goodstart","sortfactorsg")
  `[<-`(x,grep("## Optional parameter lines",x)+1,
        paste("    ",optstring,sep=""))
}
repLims <- function(x,i) {
  if( i > 1 ) {
    newval <- 0.5
    matched <- grep("##  \"lims\"",x)+1
    `[<-`(x,matched,
          sub("10\\.00000",paste(" ",format(newval,nsmall=5),sep=""),
              x[matched]))
  } else {
    return(x)
  }
}
## combined function
runNum <- function(x,i) repLims(repOpts(x,i),i)

## +++++ run PMF +++++
execPMF <- function(i,x) {
  ## INI file
  inifile <- runNum(with(x[i,],jointrepl(nFactors,FPEAK,Seed)),i)
  cat(inifile,sep="\n",file="MYPMF.INI")
  ## run PMF
  output <- system("c:/PMF/pmf2wopt.exe MYPMF",intern=TRUE)
  cat(output,sep="\n",file="")
  ## save results
  ## I do not save MISC.TXT since it is huge
  nms <- rownames(x)[i]
  try(dir.create(nms),TRUE)
  invisible(sapply(c("G_FACTOR.TXT","F_FACTOR.TXT"),
                   function(x) file.copy(x,file.path(nms,x),
                                         overwrite=TRUE)))
  cat(output,sep="\n",file=file.path(nms,"output.txt"))
  return(NULL)
}
## +++++

