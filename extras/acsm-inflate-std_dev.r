origpath <- "runs/ACSM10"
newpath <- "runs/ACSM12"
scaleby <- sqrt(3.6)

dir.create(newpath)
for( x in c("matrix.dat","samples.txt","variables.txt") ) 
  file.copy(file.path(origpath,x),newpath)
S <- as.matrix(read.table(file.path(origpath,"std_dev.dat"),colClasses="numeric"))
write.table(formatC(S*scaleby,digits=5,format="g"),file.path(newpath,"std_dev.dat"),
            sep="\t",quote=FALSE,row=FALSE,col=FALSE)
