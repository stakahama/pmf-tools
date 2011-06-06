### --- this is to reduce size of output pdf but compression is recommended ---
### --- this section can be deleted to call the original image and image.plot ---

image <- function (x = seq(0, 1, length.out = nrow(z)),
                   y = seq(0, 1, length.out = ncol(z)), z,
                   zlim = NULL,#range(z[is.finite(z)]), 
                   xlim = range(x), ylim = range(y), col = heat.colors(12), 
                   add = FALSE, xaxs = "i", yaxs = "i", xlab, ylab, breaks, 
                   oldstyle = FALSE, ...,interpolate=TRUE) 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        if (is.null(dim(x))) 
          stop("argument must be matrix-like")
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
      if (missing(xlab)) 
        xlab <- ""
      if (missing(ylab)) 
        ylab <- ""
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    xn <- deparse(substitute(x))
    if (missing(xlab)) 
      xlab <- paste(xn, "x", sep = "$")
    if (missing(ylab)) 
      ylab <- paste(xn, "y", sep = "$")
    y <- x$y
    x <- x$x
  }
  else {
    if (missing(xlab)) 
      xlab <- if (missing(x)) 
        ""
      else deparse(substitute(x))
    if (missing(ylab)) 
      ylab <- if (missing(y)) 
        ""
      else deparse(substitute(y))
  }
  if (any(!is.finite(x)) || any(!is.finite(y))) 
    stop("'x' and 'y' values must be finite and non-missing")
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  if (!is.matrix(z)) 
    stop("'z' must be a matrix")
  if (length(x) > 1 && length(x) == nrow(z)) {
    dx <- 0.5 * diff(x)
    x <- c(x[1L] - dx[1L], x[-length(x)] + dx, x[length(x)] + 
           dx[length(x) - 1])
  }
  if (length(y) > 1 && length(y) == ncol(z)) {
    dy <- 0.5 * diff(y)
    y <- c(y[1L] - dy[1L], y[-length(y)] + dy, y[length(y)] + 
           dy[length(y) - 1])
  }
  if (!add) 
    plot(NA, NA, xlim = xlim, ylim = ylim, type = "n", xaxs = xaxs, 
         yaxs = yaxs, xlab = xlab, ylab = ylab, ...)
  if (length(x) <= 1) 
    x <- par("usr")[1L:2]
  if (length(y) <= 1) 
    y <- par("usr")[3:4]
  if (length(x) != nrow(z) + 1 || length(y) != ncol(z) + 1) 
    stop("dimensions of z are not length(x)(-1) times length(y)(-1)")
  ## st
  nc <- length(col)
  if(is.null(zlim)) {
    zi <- as.raster(structure(col[cut(z,nc)],dim=dim(z)))
  } else {
    if (!missing(zlim) && (any(!is.finite(zlim)) || diff(zlim) < 
                           0))
      stop("invalid z limits")
    if (diff(zlim) == 0)
      zlim <- if (zlim[1L] == 0) 
        c(-1, 1)
      else zlim[1L] + c(-0.4, 0.4) * abs(zlim[1L])
    if(missing(breaks))
      breaks <- seq(zlim[1],zlim[2],,nc+1)
    if (any(!is.finite(breaks)))
      stop("breaks must all be finite")
    zi <- as.raster(structure(col[factor(findInterval(z,breaks),levels=1:nc)],
                                dim=dim(z)))
  }
  rasterImage(zi[seq(nrow(zi),1,-1),],
              xlim[1],ylim[1],xlim[2],ylim[2],
              interpolate=interpolate)
}

## image <- function(x,y=NULL,z=NULL,col=heat.colors(12),interpolate=TRUE,...) {
##   force(x); force(y); force(z); force(col)
##   args <- list(...)
##   if( is.null(y) && is.null(z) )
##     if( class(x) %in% "list" ) {
##       if( all(names(x) %in% letters[24:26]) ) {
##         z <- x$z
##         y <- x$y
##         x <- x$x
##       } else if( length(x)==3 ) {
##         z <- x[[3]]
##         y <- x[[2]]
##         x <- x[[1]]
##       } else {
##         print("check structure")
##       }
##     } else if( class(x)=="matrix" ) {
##       z <- x
##       y <- 1:nrow(z)
##       x <- 1:ncol(z)
##     }
##   z <- base::t(z)
##   n <- length(col)
##   dx <- diff(range(x))/(length(x)+1)
##   dy <- diff(range(y))/(length(y)+1)
##   xlim <- c(head(x,1),tail(x,1))+dx*c(-1,1)
##   ylim <- c(head(y,1),tail(y,1))+dy*c(-1,1)
##   if(is.null(args$add) || !args$add) {
##     plot.new()
##     do.call(plot.window,
##             list(xlim=xlim,
##                  ylim=ylim,
##                  asp=if(is.null(args$asp)) NA else args$asp,
##                  xaxs="i",yaxs="i"))
##   }
##   zlim <- args$zlim
##   if(is.null(zlim)) {
##     zcol <- as.raster(structure(col[cut(z,n)],dim=dim(z)))
##   } else {
##     zcol <- as.raster(structure(col[factor(findInterval(z,seq(zlim[1],zlim[2],,n+1)),
##                                            levels=1:n)],
##                                 dim=dim(z)))
##   }
##   rasterImage(zcol[seq(nrow(zcol),1,-1),],
##               xlim[1],ylim[1],xlim[2],ylim[2],
##               interpolate=interpolate)
## }

attach(NULL,name="imageplot",pos=2)
assign("image.plot",fields::image.plot,as.environment("imageplot"))
evalq(environment(image.plot) <- as.environment("imageplot"),
      as.environment("imageplot"))
