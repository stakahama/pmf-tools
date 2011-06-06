source("functions/classify.r")
source("userinputs.r")
source("ftir-inputs.r")


sc <- function(x,y) {
  i <- upper.tri(x)
  lab <- outer(rownames(x),colnames(x),paste,sep="-")
  mat <- `rownames<-`(cbind(x[i],y[i]),lab[i])
  mat[order(mat[,2],decreasing=TRUE),]
}

runno <- 26

F <- as.matrix(read.table(file.path(runnum(runno),"F_FACTOR.TXT")))
G <- as.matrix(read.table(file.path(runnum(runno),"G_FACTOR.TXT")))

sc(cor(G),cor(t(F)))

