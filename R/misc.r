rimage <- function(x = seq(0, 1, length.out = nrow(z)),
                   y = seq(0, 1, length.out = ncol(z)),
                   z, zlim=NULL, col=grey(0:255/255),
                   low="blue", up="gold", NAcolor="red",...){
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
        y <- seq.int(0, 1, length.out = ncol(z))
      }
    }
  }
   if(is.null(zlim)){
     zlim<-range(z)
   } else {
     z[z<zlim[1]] <- zlim[1]
     z[z>zlim[2]] <- zlim[2]
     col <- c(low,col,up)
   }
   image(x, y, z, zlim=zlim, col=col, ...)
   if(any(is.na(z))) image(x,y,is.na(z),col=c(NA,NAcolor),add=TRUE)
}

hg1f1 <- function(a, b, z){
  ##
  ##  Confluent Hypergeometric 1F1 (a,b scalar, z vector)
  ##  rel accuracy 1e-13 for z in -1400:700 for a=-.5, .5
  ##  rel accuracy 2e-4 for z < -1400 for a=-.5, .5
  ##
  n <- length(z)
  z[is.na(z)] <- -1e20
  z[is.infinite(z)] <- 1e-20
  .Fortran(C_hg1f1,
           as.double(a),
           as.double(b),
           as.double(z),
           as.integer(n),
           fz = double(n))$fz
}

getnlspars <- function (object) {
  r <- as.vector(object$m$resid())
  w <- object$weights
  n <- if (!is.null(w))
    sum(w > 0)
  else length(r)
  param <- coef(object)
  pnames <- names(param)
  p <- length(param)
  rdf <- n - p
  resvar <- if (rdf <= 0)
    NaN
  else deviance(object)/rdf
  Rmat <- object$m$Rmat()
  XtX <- t(Rmat)%*%Rmat
  dimnames(XtX) <- list(pnames, pnames)
  ans <- list(formula = formula(object), residuals = r, sigma = sqrt(resvar),
              df = c(p, rdf), XtX = XtX, call = object$call,
              convInfo = object$convInfo, control = object$control,
              na.action = object$na.action, coefficients = param)
  ans
}
