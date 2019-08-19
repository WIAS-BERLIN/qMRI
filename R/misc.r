extract.default <- function(x, ...) awsMethods::extract(x)
extract.MPMData <- function(x,what, ...){
  onames <- names(x)
  select <- what[what %in% onames]
  # return single component or list with selected components
    z <- x[select]
    sdim <- x$sdim
    mask <- x$mask
    if(any(dim(x$ddata)[-1]!=sdim)){
#      data are reduced to only contain voxel within mask, need to expand
       if("ddata" %in% select){
          ddata <- array(0,c(x$nFiles, prod(sdim)))
          ddata[, mask] <- x$ddata
          dim(ddata) <- c(x$nFiles, sdim)
          z[["ddata"]] <- ddata
       }
    }
  invisible(if(length(select)==1) z[[select]] else z[select])
}

extract.ESTATICSModel <- function(x,what, ...){
  onames <- names(x)
  select <- what[what %in% onames]
  z <- x[select]
  sdim <- x$sdim
  mask <- x$mask
  if(any(dim(x$modelCoeff)[-1]!=sdim)){
#      data are reduced to only contain voxel within mask, need to expand
     if("modelCoeff" %in% select){
        modelCoeff <- array(0,c(x$model+2, prod(sdim)))
        modelCoeff[, mask] <- x$modelCoeff
        dim(modelCoeff) <- c(x$model+2, sdim)
        z[["modelCoeff"]] <- modelCoeff
     }
     if("invCov" %in% select){
        invCov <- array(0,c(x$model+2, x$model+2, prod(sdim)))
        invCov[, , mask] <- x$invCov
        dim(invCov) <- c(x$model+2, x$model+2, sdim)
        z[["invCov"]] <- invCov
     }
     if("isThresh" %in% select){
        isThresh <- array(0,prod(sdim))
        isThresh[mask] <- x$isThresh
        dim(isThresh) <- sdim
        z[["isThresh"]] <- isThresh
     }
     if("isConv" %in% select){
        isConv <- array(0,prod(sdim))
        isConv[mask] <- x$isConv
        dim(isConv) <- sdim
        z[["isConv"]] <- isConv
     }
     if("sigma" %in% select){
        sigma <- array(0,prod(sdim))
        sigma[mask] <- x$sigma
        dim(sigma) <- sdim
        z[["sigma"]] <- sigma
     }
     if("rsigma" %in% select){
        sigma <- array(0,prod(sdim))
        sigma[mask] <- x$rsigma
        dim(sigma) <- sdim
        z[["rsigma"]] <- sigma
     }
  }
  invisible(if(length(select)==1) z[[select]] else z[select])
}

extract.sESTATICSModel <- function(x,what, ...){
  onames <- names(x)
  select <- what[what %in% onames]
  z <- x[select]
  sdim <- x$sdim
  mask <- x$mask
  if(any(dim(x$modelCoeff)[-1]!=sdim)){
#      data are reduced to only contain voxel within mask, need to expand
     if("modelCoeff" %in% select){
        modelCoeff <- array(0,c(x$model+2, prod(sdim)))
        modelCoeff[, mask] <- x$modelCoeff
        dim(modelCoeff) <- c(x$model+2, sdim)
        z[["modelCoeff"]] <- modelCoeff
     }
     if("invCov" %in% select){
        invCov <- array(0,c(x$model+2, x$model+2, prod(sdim)))
        invCov[, , mask] <- x$invCov
        dim(invCov) <- c(x$model+2, x$model+2, sdim)
        z[["invCov"]] <- invCov
     }
     if("isThresh" %in% select){
        isThresh <- array(0,prod(sdim))
        isThresh[mask] <- x$isThresh
        dim(isThresh) <- sdim
        z[["isThresh"]] <- isThresh
     }
     if("isConv" %in% select){
        isConv <- array(0,prod(sdim))
        isConv[mask] <- x$isConv
        dim(isConv) <- sdim
        z[["isConv"]] <- isConv
     }
     if("sigma" %in% select){
        sigma <- array(0,prod(sdim))
        sigma[mask] <- x$sigma
        dim(sigma) <- sdim
        z[["sigma"]] <- sigma
     }
     if("rsigma" %in% select){
        sigma <- array(0,prod(sdim))
        sigma[mask] <- x$rsigma
        dim(sigma) <- sdim
        z[["rsigma"]] <- sigma
     }
     if("bi" %in% select){
        bi <- array(0,prod(sdim))
        bi[mask] <- x$bi
        dim(bi) <- sdim
        z[["bi"]] <- bi
     }
     if("smoothedData" %in% select){
        ddata <- array(0,c(dim(x$smoothedData)[1], prod(sdim)))
        ddata[, mask] <- x$smoothedData
        dim(ddata) <- c(dim(x$smoothedData)[1], sdim)
        z[["smoothedData"]] <- ddata
     }
  }
  invisible(if(length(select)==1) z[[select]] else z[select])
}

extract.qMaps <- function(x,what, ...){
  onames <- names(x)
  select <- what[what %in% onames]
  z <- x[select]
  sdim <- x$sdim
  mask <- x$mask
  if(length(x$R1)==sum(mask)){
#      data are reduced to only contain voxel within mask, need to expand
     if("b1Map" %in% select){
        b1Map <- array(0,prod(sdim))
        b1Map[mask] <- x$b1Map
        dim(b1Map) <- sdim
        z[["b1Map"]] <- b1Map
     }
     if("R1" %in% select){
        R1 <- array(0,prod(sdim))
        R1[mask] <- x$R1
        dim(R1) <- sdim
        z[["R1"]] <- R1
     }
     if("R2star" %in% select){
        R2star <- array(0,prod(sdim))
        R2star[mask] <- x$R2star
        dim(R2star) <- sdim
        z[["R2star"]] <- R2star
     }
     if("PD" %in% select){
        PD <- array(0,prod(sdim))
        PD[mask] <- x$PD
        dim(PD) <- sdim
        z[["PD"]] <- PD
     }
     if("MT" %in% select){
        MT <- array(0,prod(sdim))
        MT[mask] <- x$MT
        dim(MT) <- sdim
        z[["MT"]] <- MT
     }
  }
  invisible(if(length(select)==1) z[[select]] else z[select])
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
              df = c(p, rdf), XtX = XtX, invCov = XtX/resvar, call = object$call,
              convInfo = object$convInfo, control = object$control,
              na.action = object$na.action, coefficients = param)
  ans
}

setMPMmask <- function(mpmData,mask){
   if(any(dim(mask)!=mpmData$sdim)||!is.logical(mask)){
      warning("can't set new mask returning old mpmData object \n")
      return(mpmData)
   }
   ddata <- extract(mpmData,"ddata")
   dim(ddata) <- c(mpmData$nFiles,prod(mpmData$sdim))
   mpmData$ddata <- ddata[,mask]
   mpmData$mask <- mask
   mpmData$maskFile <- "none"
   mpmData
}

getnlspars2 <- function (object, sigma, ind) {
#
#   using variance estimates from data instead of RSS
#
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
  grad <- object$m$gradient()
  XtX <- t(grad)%*%grad
  sgrad <- sigma[ind] * grad
  z <- svd(sgrad)
  if(any(z$d<1e-6*max(z$d))){
     cat("singular covariance\ngradient:\n")
     print(grad)
     cat("sigma\n")
     print(sigma[ind])
  }
  z$d <- pmax(z$d,1e-6*max(z$d))
  sXtXinv <- z$v%*%diag(1/z$d^2)%*%t(z$v)
  XtXsinv <- XtX%*%sXtXinv%*%XtX
  dimnames(XtX) <- list(pnames, pnames)
  ans <- list(formula = formula(object), residuals = r, sigma = sqrt(resvar),
              df = c(p, rdf), XtX = XtX, invCov=XtXsinv,call = object$call,
              convInfo = object$convInfo, control = object$control,
              na.action = object$na.action, coefficients = param)
  ans
}
