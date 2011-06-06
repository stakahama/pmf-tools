IGORdatetime <- function (x) 
  as.chron(ISOdatetime(1904, 1, 1, 0, 0, 0, "GMT") + x)

format_chron <- function (x, fmt = "%Y-%m-%d %H:%M:%S") 
  format(as.POSIXct(x), tz = "GMT", fmt)

  append.default <- base::append
  append <- `body<-`(args(append.default),value=quote(UseMethod("append")))
  append.data.frame <-
    `body<-`(args(append.default),value=
             quote(`row.names<-`(data.frame(append.default(x,values,after),check.names=FALSE),
                                 row.names(x))))

zip <- function(...) {
  dotArgs <- list(...)
  arglist <- (if( length(dotArgs)==1 &&
                 is.list(dotArgs[[1]]) )
              dotArgs[[1]] else dotArgs)
  do.call(function(...) Map(list,...),arglist)
}
