extract.default <- function(x, ...) awsMethods::extract(x)
extract.MPMData <- function(x,what, ...){
  onames <- names(x)
  select <- what[what %in% onames]
  # return single component or list with selected components
  invisible(if(length(select)==1) x[[select]] else x[select])
}

extract.ESTATICSModel <- function(x,what, ...){
  onames <- names(x)
  select <- what[what %in% onames]
  invisible(if(length(select)==1) x[[select]] else x[select])
}

extract.qMaps <- function(x,what, ...){
  onames <- names(x)
  select <- what[what %in% onames]
  invisible(if(length(select)==1) x[[select]] else x[select])
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
