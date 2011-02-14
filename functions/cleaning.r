####################
## PMF execution and postprocessing program
## ~cleaning.r~
## $Rev$
## Feb. 2011
## Satoshi Takahama (stakahama@ucsd.edu)
####################


cleanbl <- function(w,a) {
  negp <- max(w[a < 0 & w < 3000 & w > 2000])
  if( is.finite(negp) ) a[w <= negp & w > 2000] <- 0
  negp <- min(w[a < 0 & w > 3500])
  if( is.finite(negp) ) a[w >= negp] <- 0
  data.frame(wavenumber=w,baselined=a)
}
