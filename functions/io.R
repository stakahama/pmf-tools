

library(RJSONIO)

read.args <- function(filename) {
  args <- fromJSON(filename)
  args[["FPEAK"]] <- do.call(seq,as.list(args[["FPEAK"]][c(1,3,2)]))
  if(is.null(args[["FOLDER"]]))
    args[["FOLDER"]] <- path.expand(dirname(filename))
  args
}
