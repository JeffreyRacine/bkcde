## Kernel functions for bkcde package

## pdf.kernel.bk() is the doubly truncated Gaussian boundary kernel function from
## Racine et al 2024

pdf.kernel.bk <- function(x,X,h,a=-Inf,b=Inf,denom=NULL) {
  ## Checking for bounds involves a bit of overhead (20%), so here we presume a
  ## check is performed outside of this function - make sure this is the case!
  z <- (x-X)/h
  if(is.null(denom)) {
    p1 <- if(is.infinite(b)) 1 else pnorm((b-x)/h)
    p2 <- if(is.infinite(a)) 0 else pnorm((a-x)/h)
    denom <- h*(p1-p2)
  }
  (exp(-0.5 * z*z) * 0.3989422804014326779) / denom
}

## cdf.kernel.bk() is the doubly truncated Gaussian boundary kernel function for
## the distribution function from Racine et al 2024

cdf.kernel.bk <- function(x,X,h,a=-Inf,b=Inf,denom=NULL,pnorm.a=NULL) {
  ## Checking for bounds involves a bit of overhead (20%), so here we presume a
  ## check is performed outside of this function - make sure this is the case!
  if(is.null(pnorm.a)) pnorm.a <- if(is.infinite(a)) 0 else pnorm((a-X)/h)
  if(is.null(denom)) {
    pnorm.b <- if(is.infinite(b)) 1 else pnorm((b-X)/h)
    denom <- pnorm.b - pnorm.a
  }
  (pnorm((x-X)/h)-pnorm.a)/denom
}
