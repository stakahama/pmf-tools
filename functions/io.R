

library(RJSONIO)

readArgs <- function(filename) {
  args <- fromJSON(filename)
  args[["FPEAK"]] <- do.call(seq,as.list(args[["FPEAK"]][c(1,3,2)]))
  if(is.null(args[["FOLDER"]]))
    args[["FOLDER"]] <- path.expand(dirname(filename))
  args
}


convertFactor <- function(df, columns) {
  for(i in columns) {
    x <- df[[i]]
    df[[i]] <- as.integer(levels(x))[x]
  }
  df
}
