
datapath <- "data"
mat <- as.matrix(read.table(file.path(datapath,"matrix.dat")))
std <- as.matrix(read.table(file.path(datapath,"std_dev.dat")))
wave <- scan(file.path(datapath,"wavenumbers.txt"),quiet=TRUE)
samp <- readLines(file.path(datapath,"samples.txt"))


folder <- "runs/TJ3"
dir.create(folder)
std[] <- 1/abs(mat)

colsplit <- function(x) split(x,col(x))
std[] <- mapply(function(x,y) if(any(is.finite(y))) replace(y,!is.finite(y),max(x)*4) else y,
                colsplit(mat), colsplit(std))

ix <- apply(mat,1,function(x,w)
            all(x[w > 3000 & w < 3400] > 0),
            wave)
## jx <- apply(std,2,function(x) all(x > 0))
jx <- apply(std,2,function(x) all(is.finite(x)))

write.table <- function(...)
  utils::write.table(...,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(format(mat[ix,jx],nsmall=5),file.path(folder,"matrix.dat"))
write.table(format(std[ix,jx],nsmall=5),file.path(folder,"std_dev.dat"))
writeLines(format(wave[jx],nsmall=3),file.path(folder,"wavenumbers.txt"))
writeLines(samp[ix],file.path(folder,"samples.txt"))

## matplot(wave[jx],t(mat[ix,jx]),type="l",lty=1,col=8)

matplot(wave[jx],t(std[ix,jx]),type="l",lty=1,col=8)

