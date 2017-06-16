readMPMData  <-  function(t1Files  = NULL,
                          pdFiles  = NULL,
                          mtFiles  = NULL,
                          maskFile = NULL,
                          sdim     = NULL,
                          TR       = NULL,
                          TE       = NULL,
                          FA       = NULL,
                          verbose  = TRUE) {

  ## we need at least T1w and PDw files
  if (is.null(t1Files)) stop("vector of T1 files required")
  if (is.null(pdFiles)) stop("vector of PD files required")
  ## TODO: test whether there are enough files for the model?

  if (is.null(sdim)) stop("need spatial dimensionality of the data")
  if (!is.numeric(sdim) | length(sdim) != 3) stop("need exactly three numbers for spatial dimensions")

  ## select the model according to the existence of MTw files
  model <- if (is.null(mtFiles)) {
    1L # the simple model without MT inclusion
  } else {
    2L # the model including MT
  }

  ## count the number of data volumes (if (is.null(mtFiles)) length(mtFiles) == 0)
  nFiles <- length(t1Files) + length(mtFiles) + length(pdFiles)

  # the array for the data itself
  ddata <- array(0, c(nFiles, prod(sdim)))

  ## for each files we have a TR, TE, and flip angle (FA)
  if (is.null(TR) || is.null(TE) || is.null(FA)) {
    TR <- TE <- FA <- numeric(nFiles)
    readParameterFlag <- TRUE
  } else {
    if (length(TR) != nFiles) stop("not enough TR value, need as many as file!")
    if (length(TE) != nFiles) stop("not enough TE value, need as many as file!")
    if (length(FA) != nFiles) stop("not enough FA value, need as many as file!")
    readParameterFlag <- FALSE
  }

  ## ... now we read all data volumes and extract the TR, TE, and FA values for each ...
  ii <- 1
  ## ... for all T1 volumes ...
  if (verbose) cat("reading T1 files\n")
  if (verbose) pb <- txtProgressBar(min = 0, max = length(t1Files), style = 3)
  for (i in 1:length(t1Files)) {
    ds <- readNIfTI(t1Files[i], reorient = FALSE)
    ddata[ii, ] <- ds
    if(readParameterFlag) {
      ## IMPORTANT: This is special to Siawoosh data
      res <- str_match_all(ds@descrip, "([[:alpha:]]{2})=([.0123456789]+)([[:alpha:]]{2,})")
      for (nn in 1:dim(res[[1]])[1]) {
        if (res[[1]][nn, 2] == "TR") TR[ii] <- as.numeric(res[[1]][nn, 3])
        if (res[[1]][nn, 2] == "TE") TE[ii] <- as.numeric(res[[1]][nn, 3])
        if (res[[1]][nn, 2] == "FA") FA[ii] <- as.numeric(res[[1]][nn, 3])
      }
    }
    ii <- ii + 1
    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) close(pb)
  ## ... for all MT volumes ...
  if (model == 2) {
    if (verbose) cat("reading MT files\n")
    if (verbose) pb <- txtProgressBar(min = 0, max = length(mtFiles), style = 3)
    for (i in 1:length(mtFiles)) {
      ds <- readNIfTI(mtFiles[i], reorient = FALSE)
      ddata[ii, ] <- ds
      if(readParameterFlag) {
        ## IMPORTANT: This is special to Siawoosh data
        res <- str_match_all(ds@descrip, "([[:alpha:]]{2})=([.0123456789]+)([[:alpha:]]{2,})")
        for (nn in 1:dim(res[[1]])[1]) {
          if (res[[1]][nn, 2] == "TR") TR[ii] <- as.numeric(res[[1]][nn, 3])
          if (res[[1]][nn, 2] == "TE") TE[ii] <- as.numeric(res[[1]][nn, 3])
          if (res[[1]][nn, 2] == "FA") FA[ii] <- as.numeric(res[[1]][nn, 3])
        }
      }
      ii <- ii + 1
      if (verbose) setTxtProgressBar(pb, i)
    }
  }
  if (verbose) close(pb)
  ## .. and for all PD volumes ...
  if (verbose) cat("reading PD files\n")
  if (verbose) pb <- txtProgressBar(min = 0, max = length(pdFiles), style = 3)
  for (i in 1:length(pdFiles)) {
    ds <- readNIfTI(pdFiles[i], reorient = FALSE)
    ddata[ii, ] <- ds
    if(readParameterFlag) {
      ## IMPORTANT: This is special to Siawoosh data
      res <- str_match_all(ds@descrip, "([[:alpha:]]{2})=([.0123456789]+)([[:alpha:]]{2,})")
      for (nn in 1:dim(res[[1]])[1]) {
        if (res[[1]][nn, 2] == "TR") TR[ii] <- as.numeric(res[[1]][nn, 3])
        if (res[[1]][nn, 2] == "TE") TE[ii] <- as.numeric(res[[1]][nn, 3])
        if (res[[1]][nn, 2] == "FA") FA[ii] <- as.numeric(res[[1]][nn, 3])
      }
    }
    ii <- ii + 1
    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) close(pb)
  dim(ddata) <- c(nFiles, sdim)
  ## ... done!

  # the mask array
  ## TODO: set mask TRUE also if maskFile does not exist or cannot be read by readNIFTI()
  if (is.null(maskFile)) {
    if (verbose) cat("found no mask file, setting mask TRUE everywhere")
    mask <- array(TRUE, sdim)
  } else {
    if (verbose) cat("reading mask file ... ")
    mask <- as.logical(readNIfTI(maskFile, reorient = FALSE))
    dim(mask) <- sdim
    if (verbose) cat("done\n")
  }

  invisible(list(ddata = ddata,
                 sdim = sdim,
                 nFiles = nFiles,
                 t1Files = t1Files,
                 pdFiles = pdFiles,
                 mtFiles = mtFiles,
                 model = model,
                 maskFile = maskFile,
                 mask = mask,
                 TR = TR,
                 TE = TE,
                 FA = FA))
}

estimateSigma <- function(magnitude,phase,mask,kstar=20,kmin=8,hsig=5,lambda=12,verbose=TRUE){
  ## kmin = 10 corresponds to an initial bandwidth of 1.47 giving positive weight to direct neighbors and 
  ## 2D diagonal neigbors 
  args <- sys.call(-1)
  sdim <- dim(mask)
  if(!is.numeric(magnitude)){
     if (verbose) cat("reading Magnitude file ... ")
     R <- readNIfTI(magnitude, reorient = FALSE)
  } else {
  R <- magnitude
  }
  if(!is.numeric(phase)){
    if (verbose) cat("reading Phase file ... ")
     Ph <- readNIfTI(phase, reorient = FALSE)
  } else {
     Ph <- phase
  }
  ComplImg <- array(0,c(2,sdim))
  ComplImg[1,,,] <- R*cos(Ph)
  ComplImg[2,,,] <- R*sin(Ph)
  ## find the number of usable cores
  mc.cores <- setCores(, reprt = FALSE)
##
##  start smoothing and variance estimation
##
  n <- prod(sdim)
  lambda0 <- 1e40
  sigma2 <- array(1e10,sdim)
# just inilitialize with something large, first step is nonadaptive due to lambda0 
  k <- kmin
  hmax <- 1.25^(kstar/3)
  ## preparations for median smoothing
  nwmd <- (2*as.integer(hsig)+1)^3
  parammd <- .Fortran("paramw3",
                      as.double(hsig),
                      as.double(c(1,1)),
                      ind=integer(3*nwmd),
                      w=double(nwmd),
                      n=as.integer(nwmd),
                      PACKAGE = "qMRI")[c("ind","w","n")]
  nwmd <- parammd$n
  parammd$ind <- parammd$ind[1:(3*nwmd)]
  dim(parammd$ind) <- c(3,nwmd)
  
  if (verbose) pb <- txtProgressBar(min = 0, max = kstar-kmin+1, style = 3)
  bi <- array(1,sdim)
  zobj <- list(theta=ComplImg, bi=bi)
  if (verbose) {
    mae <- NULL
    protocol <- matrix("", kstar-kmin+1, 1, dimnames = list(paste("step", kmin:kstar), "protocol"))
  }
  while (k <= kstar) {
    ## determine the actual bandwidth for this step
    hakt <- gethani(1, 1.25*hmax, 2, 1.25^k, c(1,1), 1e-4)
    
    ## we need the (approx.) size of the weigthing scheme array
    dlw <- (2*trunc(hakt/c(1, 1, 1))+1)[1:3]
    
    ## perform the actual adaptive smoothing
    zobj <- .Fortran("vaws2",
                     as.double(ComplImg),
                     as.logical(mask),
                     as.integer(2),
                     as.integer(sdim[1]),
                     as.integer(sdim[2]),
                     as.integer(sdim[3]),
                     hakt = as.double(hakt),
                     as.double(lambda0),
                     as.double(zobj$theta),
                     as.double(sigma2),
                     bi = as.double(zobj$bi),
                     theta = double(2*n),
                     sigma2 = double(n),
                     as.integer(mc.cores),
                     double(prod(dlw)),
                     as.double(c(1,1)),
                     double(2 * mc.cores),
                     PACKAGE = "qMRI")[c("bi", "theta", "hakt","sigma2")]
    ##
    ##  now get local median variance estimates
    ##
    dim(zobj$sigma2) <- sdim
    sigma2 <- .Fortran("mediansm",
                      as.double(zobj$sigma2),
                      as.logical(mask),
                      as.integer(sdim[1]),
                      as.integer(sdim[2]),
                      as.integer(sdim[3]),
                      as.integer(parammd$ind),
                      as.integer(nwmd),
                      double(nwmd*mc.cores), # work(nw,nthreds)
                      as.integer(mc.cores),
                      sigma2n = double(n),
                      PACKAGE = "qMRI")$sigma2n/0.6931
# sigma2n containes sum of 2 independent squared residuals
# 0.6931 approximates  median \chi_2 /2  
# needed to get correct results 
    ## use maximum ni
    bi <- zobj$bi <- pmax(bi, zobj$bi)
    
    ## some verbose stuff
    if (verbose) {
      protocol[k-kmin+1,1] <- paste("bandwidth: ", signif(hakt, 3),
                           "sigma: mean: ", signif(sqrt(mean(sigma2[mask])),3),
                           "median: ", signif(sqrt(median(sigma2[mask])),3),
                           "sd: ", signif(sd(sqrt(sigma2[mask])),3),
                           "median(bi):", signif(median(zobj$bi[mask]),3),
                           "max(bi):", signif(max(zobj$bi[mask]),3))
      setTxtProgressBar(pb, k-kmin+1)
    }
    
    ## go for next iteration
    k <- k+1
    lambda0 <- lambda
    gc()
  }
  dim(zobj$theta) <- c(2,sdim)
  # return estimated parameters of rician distribution
  z <-  list(sigma=array(sqrt(sigma2),sdim),
       theta=array(sqrt(zobj$theta[1,,,]^2+zobj$theta[2,,,]^2),sdim),
       sigmal=array(sqrt(zobj$sigma2),sdim),mask=mask, 
       protocol=protocol,args=args)
  class(z) <- "sigmaEstSENSE"
  z
}

medianFilterSigma <- function(obj,hsig=10,mask=NULL){
   if(class(obj)=="sigmaEstSENSE"){
      sigma2 <- obj$sigmal^2
      mask <- obj$mask
   } else {
      sigma2 <- obj^2
   }
   sdim <- dim(sigma2)
   n <- prod(sdim)
   if(length(sdim)!=3) stop("obj needs to be of class 'array' (3D) or 'sigmaEstSENSE'")
   if(is.null(mask)) mask <- array(TRUE,sdim)
   if(any(dim(mask)!=sdim)) stop("dimensions do not coinside")
   nwmd <- (2*as.integer(hsig)+1)^3
   parammd <- .Fortran("paramw3",
                       as.double(hsig),
                       as.double(c(1,1)),
                       ind=integer(3*nwmd),
                       w=double(nwmd),
                       n=as.integer(nwmd),
                       PACKAGE = "qMRI")[c("ind","w","n")]
   nwmd <- parammd$n
   parammd$ind <- parammd$ind[1:(3*nwmd)]
   dim(parammd$ind) <- c(3,nwmd)
   mc.cores <- setCores(, reprt = FALSE)
   sigma2 <- .Fortran("mediansm",
                      as.double(sigma2),
                      as.logical(mask),
                      as.integer(sdim[1]),
                      as.integer(sdim[2]),
                      as.integer(sdim[3]),
                      as.integer(parammd$ind),
                      as.integer(nwmd),
                      double(nwmd*mc.cores), # work(nw,nthreds)
                      as.integer(mc.cores),
                      sigma2n = double(n),
                      PACKAGE = "qMRI")$sigma2n/0.6931
   dim(sigma2) <- sdim
   if(class(obj)=="sigmaEstSENSE"){
      obj$sigma <- sqrt(sigma2)
      obj$hsig <- hsig
   } else {
      obj <- sqrt(sigma2)
   }
   obj
}
# estimateESTATICS <- function(mpmdata,
#                              TEScale = 100,
#                              dataScale = 1000,
#                              verbose = TRUE) {
# 
#   ## this is our design ...
#   if (mpmdata$model == 2) {
#     xmat <- matrix(0, mpmdata$nFiles, 4)
#     xmat[1:length(mpmdata$t1Files), 1] <- 1
#     xmat[(length(mpmdata$t1Files)+1):(length(mpmdata$t1Files)+length(mpmdata$mtFiles)), 2] <- 1
#     xmat[(length(mpmdata$t1Files)+length(mpmdata$mtFiles)+1):mpmdata$nFiles, 3] <- 1
#     xmat[, 4] <- mpmdata$TE / TEScale
#     ## ... for our model in qflashpl() ...
#     ## S_{T1} = par[1] * exp(- par[4] * TE)
#     ## S_{MT} = par[2] * exp(- par[4] * TE)
#     ## S_{PD} = par[3] * exp(- par[4] * TE)
#   } else {
#     xmat <- matrix(0, mpmdata$nFiles, 3)
#     xmat[1:length(mpmdata$t1Files), 1] <- 1
#     xmat[(length(mpmdata$t1Files)+1):mpmdata$nFiles, 2] <- 1
#     xmat[, 3] <- mpmdata$TE / TEScale
#     ## ... for our model in qflashpl2() ...
#     ## S_{T1} = par[1] * exp(- par[3] * TE)
#     ## S_{PD} = par[2] * exp(- par[3] * TE)
#   }
#   if (verbose) {
#     cat("Design of the model:\n")
#     print(xmat)
#   }
# 
#   ## starting value for R* estimate
#   R2star <- 0.05 * TEScale
#   indT1 <- order(mpmdata$TE[as.logical(xmat[, 1])])[1]
#   if (mpmdata$model == 2) {
#     indMT <- order(mpmdata$TE[as.logical(xmat[, 2])])[1] + sum(xmat[, 1])
#     indPD <- order(mpmdata$TE[as.logical(xmat[, 3])])[1] + sum(xmat[, 1]) + sum(xmat[, 2])
#     npar <- 4
#   } else {
#     indPD <- order(mpmdata$TE[as.logical(xmat[, 2])])[1] + sum(xmat[, 1])
#     npar <- 3
#   }
# 
#   isConv <- array(FALSE, mpmdata$sdim)
#   modelCoeff <- array(0, c(npar, mpmdata$sdim))
#   invCov <- array(0, c(npar, npar, mpmdata$sdim))
# 
#   ## now perform the voxelwise regression
#   if (verbose) Sys.time()
#   for (z in 1:mpmdata$sdim[3]){
#     for (y in 1:mpmdata$sdim[2]) {
#       for (x in 1:mpmdata$sdim[1]) {
#         if (mpmdata$mask[x, y, z]) {
#           ivec  <- mpmdata$ddata[, x, y, z] / dataScale
#           if (mpmdata$model == 2) {
#             ## full ESTATICS model
#             th <- c(ivec[indT1] * exp(-xmat[indT1, 4] * R2star), # par[1]
#                     ivec[indMT] * exp(-xmat[indMT, 4] * R2star), # par[2]
#                     ivec[indPD] * exp(-xmat[indPD, 4] * R2star), # par[3]
#                     R2star)                                      # par[4]
#             res <- try(nls(ivec ~ qflashpl(par, xmat),
#                            start = list(par = th),
#                            control = list(maxiter = 200,
#                                           warnOnly = TRUE)))
#             if(class(res) == "try-error" || !res$convInfo$isConv || any(coefficients(res) < 0))
#               res <- nls(ivec ~ qflashpl(par, xmat),
#                          start = list(par = th),
#                          algorithm = "port",
#                          control = list(warnOnly = TRUE,
#                                         printEval = TRUE),
#                          lower = rep(0, 4))
#           } else {
#             ## reduced ESTATICS model without MT
#             th <- c(ivec[indT1] * exp(-xmat[indT1, 3] * R2star), # par[1]
#                     ivec[indPD] * exp(-xmat[indPD, 3] * R2star), # par[2]
#                     R2star)                                      # par[3]
#             res <- try(nls(ivec ~ qflashpl2(par, xmat),
#                            start = list(par = th),
#                            control = list(maxiter = 200,
#                                           warnOnly = TRUE)))
#             if(class(res) == "try-error" || !res$convInfo$isConv || any(coefficients(res) < 0))
#               res <- nls(ivec ~ qflashpl2(par, xmat),
#                          start = list(par = th),
#                          algorithm = "port",
#                          control = list(warnOnly = TRUE,
#                                         printEval = TRUE),
#                          lower = rep(0, 3))
#           }
#           sres <- summary(res)
#           isConv[x, y, z] <- res$convInfo$isConv
#           modelCoeff[, x, y, z] <- sres$coefficients[, 1]
#           if (sres$sigma != 0) {
#             invCov[, , x, y, z] <- solve(sres$cov.unscaled) / sres$sigma^2
#           }
#         }
#       }
#     }
#     if (verbose) cat(z, format(Sys.time()), "\n")
#   }
#   if (verbose) Sys.time()
# 
#   invisible(list(modelCoeff = modelCoeff,
#                  invCov = invCov,
#                  isConv = isConv,
#                  sdim = mpmdata$sdim,
#                  nFiles = mpmdata$nFiles,
#                  t1Files = mpmdata$t1Files,
#                  pdFiles = mpmdata$pdFiles,
#                  mtFiles = mpmdata$mtFiles,
#                  model = mpmdata$model,
#                  maskFile = mpmdata$maskFile,
#                  mask = mpmdata$mask,
#                  TR = mpmdata$TR,
#                  TE = mpmdata$TE,
#                  FA = mpmdata$FA,
#                  TEScale = TEScale,
#                  dataScale = dataScale))
# 
#   ## END function estimateESTATICS()
# }

estimateESTATICS <- function(mpmdata,
                             TEScale = 100,
                             dataScale = 1000,
                             verbose = TRUE) {
  
  ## this is our design ...
  if (mpmdata$model == 2) {
    xmat <- matrix(0, mpmdata$nFiles, 4)
    xmat[1:length(mpmdata$t1Files), 1] <- 1
    xmat[(length(mpmdata$t1Files)+1):(length(mpmdata$t1Files)+length(mpmdata$mtFiles)), 2] <- 1
    xmat[(length(mpmdata$t1Files)+length(mpmdata$mtFiles)+1):mpmdata$nFiles, 3] <- 1
    xmat[, 4] <- mpmdata$TE / TEScale
    ## ... for our model in qflashpl() ...
    ## S_{T1} = par[1] * exp(- par[4] * TE)
    ## S_{MT} = par[2] * exp(- par[4] * TE)
    ## S_{PD} = par[3] * exp(- par[4] * TE)
    #
    #   check for voxel in mask with all zeros for a modality
    #
    zerovoxel <- as.logical(((apply(mpmdata$ddata[xmat[,1]==1,,,],2:4,sum)==0)|
                               (apply(mpmdata$ddata[xmat[,2]==1,,,],2:4,sum)==0)|
                               (apply(mpmdata$ddata[xmat[,3]==1,,,],2:4,sum)==0))*mpmdata$mask)
  } else {
    xmat <- matrix(0, mpmdata$nFiles, 3)
    xmat[1:length(mpmdata$t1Files), 1] <- 1
    xmat[(length(mpmdata$t1Files)+1):mpmdata$nFiles, 2] <- 1
    xmat[, 3] <- mpmdata$TE / TEScale
    ## ... for our model in qflashpl2() ...
    ## S_{T1} = par[1] * exp(- par[3] * TE)
    ## S_{PD} = par[2] * exp(- par[3] * TE)
    #
    #   check for voxel in mask with all zeros for a modality
    #
    zerovoxel <- as.logical(((apply(mpmdata$ddata[xmat[,1]==1,,,],2:4,sum)==0)|
                               (apply(mpmdata$ddata[xmat[,2]==1,,,],2:4,sum)==0))*mpmdata$mask)
  }
  if (verbose) {
    cat("Design of the model:\n")
    print(xmat)
  }
  #
  #   deactivate voxel in mask with all zeros for a modality
  #
  mpmdata$mask[zerovoxel] <- FALSE
  
  ## starting value for R* estimate
  R2star <- 0.05 * TEScale
  indT1 <- order(mpmdata$TE[as.logical(xmat[, 1])])[1]
  if (mpmdata$model == 2) {
    indMT <- order(mpmdata$TE[as.logical(xmat[, 2])])[1] + sum(xmat[, 1])
    indPD <- order(mpmdata$TE[as.logical(xmat[, 3])])[1] + sum(xmat[, 1]) + sum(xmat[, 2])
    npar <- 4
  } else {
    indPD <- order(mpmdata$TE[as.logical(xmat[, 2])])[1] + sum(xmat[, 1])
    npar <- 3
  }
  
  isConv <- array(FALSE, mpmdata$sdim)
  modelCoeff <- array(0, c(npar, mpmdata$sdim))
  invCov <- array(0, c(npar, npar, mpmdata$sdim))
  
  ## now perform the voxelwise regression
  if (verbose) Sys.time()
  for (z in 1:mpmdata$sdim[3]){
    for (y in 1:mpmdata$sdim[2]) {
      for (x in 1:mpmdata$sdim[1]) {
        if (mpmdata$mask[x, y, z]) {
          ivec  <- mpmdata$ddata[, x, y, z] / dataScale
          if (mpmdata$model == 2) {
            ## full ESTATICS model
            th <- c(ivec[indT1] * exp(-xmat[indT1, 4] * R2star), # par[1]
                    ivec[indMT] * exp(-xmat[indMT, 4] * R2star), # par[2]
                    ivec[indPD] * exp(-xmat[indPD, 4] * R2star), # par[3]
                    R2star)                                      # par[4]
            res <- try(nls(ivec ~ qflashpl(par, xmat),
                           start = list(par = th),
                           control = list(maxiter = 200,
                                          warnOnly = TRUE)))
            if(class(res) == "try-error" || !res$convInfo$isConv || any(coefficients(res) < 0))
              res <- nls(ivec ~ qflashpl(par, xmat),
                         start = list(par = th),
                         algorithm = "port",
                         control = list(warnOnly = TRUE,
                                        printEval = TRUE),
                         lower = rep(0, 4))
          } else {
            ## reduced ESTATICS model without MT
            th <- c(ivec[indT1] * exp(-xmat[indT1, 3] * R2star), # par[1]
                    ivec[indPD] * exp(-xmat[indPD, 3] * R2star), # par[2]
                    R2star)                                      # par[3]
            res <- try(nls(ivec ~ qflashpl2(par, xmat),
                           start = list(par = th),
                           control = list(maxiter = 200,
                                          warnOnly = TRUE)))
            if(class(res) == "try-error" || !res$convInfo$isConv || any(coefficients(res) < 0))
              res <- nls(ivec ~ qflashpl2(par, xmat),
                         start = list(par = th),
                         algorithm = "port",
                         control = list(warnOnly = TRUE,
                                        printEval = TRUE),
                         lower = rep(0, 3))
          }
          sres <- summary(res)
          isConv[x, y, z] <- res$convInfo$isConv
          modelCoeff[, x, y, z] <- sres$coefficients[, 1]
          if (sres$sigma != 0) {
            invCov[, , x, y, z] <- solve(sres$cov.unscaled) / sres$sigma^2
          }
        }
      }
    }
    if (verbose) cat(z, format(Sys.time()), "\n")
  }
  if (verbose) Sys.time()
  
  invisible(list(modelCoeff = modelCoeff,
                 invCov = invCov,
                 isConv = isConv,
                 sdim = mpmdata$sdim,
                 nFiles = mpmdata$nFiles,
                 t1Files = mpmdata$t1Files,
                 pdFiles = mpmdata$pdFiles,
                 mtFiles = mpmdata$mtFiles,
                 model = mpmdata$model,
                 maskFile = mpmdata$maskFile,
                 mask = mpmdata$mask,
                 TR = mpmdata$TR,
                 TE = mpmdata$TE,
                 FA = mpmdata$FA,
                 TEScale = TEScale,
                 dataScale = dataScale))
  
  ## END function estimateESTATICS()
}
# estimateESTATICSQL <- function(mpmdata,
#                                TEScale = 100,
#                                dataScale = 1000,
#                                sigma = NULL,
#                                L = NULL,
#                                verbose = TRUE) {
# 
#   ## this is our design ...
#   if (mpmdata$model == 2) {
#     xmat <- matrix(0, mpmdata$nFiles, 4)
#     xmat[1:length(mpmdata$t1Files), 1] <- 1
#     xmat[(length(mpmdata$t1Files)+1):(length(mpmdata$t1Files)+length(mpmdata$mtFiles)), 2] <- 1
#     xmat[(length(mpmdata$t1Files)+length(mpmdata$mtFiles)+1):mpmdata$nFiles, 3] <- 1
#     xmat[, 4] <- mpmdata$TE / TEScale
#     ## ... for our model in qflashpl() ...
#     ## S_{T1} = par[1] * exp(- par[4] * TE)
#     ## S_{MT} = par[2] * exp(- par[4] * TE)
#     ## S_{PD} = par[3] * exp(- par[4] * TE)
#   } else {
#     xmat <- matrix(0, mpmdata$nFiles, 3)
#     xmat[1:length(mpmdata$t1Files), 1] <- 1
#     xmat[(length(mpmdata$t1Files)+1):mpmdata$nFiles, 2] <- 1
#     xmat[, 3] <- mpmdata$TE / TEScale
#     ## ... for our model in qflashpl2() ...
#     ## S_{T1} = par[1] * exp(- par[3] * TE)
#     ## S_{PD} = par[2] * exp(- par[3] * TE)
#   }
#   if (verbose) {
#     cat("Design of the model:\n")
#     print(xmat)
#   }
# 
#   sigma <- sigma / dataScale
#   CL <- sigma * sqrt(pi/2) * gamma(L+0.5) / gamma(L) / gamma(1.5)
# 
#   ## starting value for R* estimate
#   R2star <- 0.05 * TEScale
#   indT1 <- order(mpmdata$TE[as.logical(xmat[, 1])])[1]
#   if (mpmdata$model == 2) {
#     indMT <- order(mpmdata$TE[as.logical(xmat[, 2])])[1] + sum(xmat[, 1])
#     indPD <- order(mpmdata$TE[as.logical(xmat[, 3])])[1] + sum(xmat[, 1]) + sum(xmat[, 2])
#     npar <- 4
#   } else {
#     indPD <- order(mpmdata$TE[as.logical(xmat[, 2])])[1] + sum(xmat[, 1])
#     npar <- 3
#   }
# 
#   isConv <- array(FALSE, mpmdata$sdim)
#   modelCoeff <- array(0, c(npar, mpmdata$sdim))
#   invCov <- array(0, c(npar, npar, mpmdata$sdim))
# 
#   ## now perform the voxelwise regression
#   if (verbose) Sys.time()
#   for (z in 1:mpmdata$sdim[3]){
#     for (y in 1:mpmdata$sdim[2]) {
#       for (x in 1:mpmdata$sdim[1]) {
#         if (mpmdata$mask[x, y, z]) {
#           ivec  <- mpmdata$ddata[, x, y, z] / dataScale
#           if (mpmdata$model == 2) {
#             ## full ESTATICS model
#             th <- c(ivec[indT1] * exp(-xmat[indT1, 4] * R2star), # par[1]
#                     ivec[indMT] * exp(-xmat[indMT, 4] * R2star), # par[2]
#                     ivec[indPD] * exp(-xmat[indPD, 4] * R2star), # par[3]
#                     R2star)                                      # par[4]
#             res <- try(nls(ivec ~ qflashplQL(par, xmat, CL, sigma, L),
#                            start = list(par = th),
#                            control = list(maxiter = 200,
#                                           warnOnly = TRUE)))
#             if(class(res) == "try-error" || !res$convInfo$isConv || any(coefficients(res) < 0))
#               res <- nls(ivec ~ qflashplQL(par, xmat, CL, sigma, L),
#                          start = list(par = th),
#                          algorithm = "port",
#                          control = list(warnOnly = TRUE,
#                                         printEval = TRUE),
#                          lower = rep(0, 4))
#           } else {
#             ## reduced ESTATICS model without MT
#             th <- c(ivec[indT1] * exp(-xmat[indT1, 3] * R2star), # par[1]
#                     ivec[indPD] * exp(-xmat[indPD, 3] * R2star), # par[2]
#                     R2star)                                      # par[3]
#             res <- try(nls(ivec ~ qflashpl2QL(par, xmat, CL, sigma, L),
#                            start = list(par = th),
#                            control = list(maxiter = 200,
#                                           warnOnly = TRUE)))
#             if(class(res) == "try-error" || !res$convInfo$isConv || any(coefficients(res) < 0))
#               res <- nls(ivec ~ qflashpl2QL(par, xmat, CL, sigma, L),
#                          start = list(par = th),
#                          algorithm = "port",
#                          control = list(warnOnly = TRUE,
#                                         printEval = TRUE),
#                          lower = rep(0, 3))
#           }
#           sres <- summary(res)
#           isConv[x, y, z] <- res$convInfo$isConv
#           modelCoeff[, x, y, z] <- sres$coefficients[, 1]
#           if (sres$sigma != 0) {
#             invCov[, , x, y, z] <- solve(sres$cov.unscaled) / sres$sigma^2
#           }
#         }
#       }
#     }
#     if (verbose) cat(z, format(Sys.time()), "\n")
#   }
#   if (verbose) Sys.time()
# 
#   invisible(list(modelCoeff = modelCoeff,
#                  invCov = invCov,
#                  isConv = isConv,
#                  sdim = mpmdata$sdim,
#                  nFiles = mpmdata$nFiles,
#                  t1Files = mpmdata$t1Files,
#                  pdFiles = mpmdata$pdFiles,
#                  mtFiles = mpmdata$mtFiles,
#                  model = mpmdata$model,
#                  maskFile = mpmdata$maskFile,
#                  mask = mpmdata$mask,
#                  TR = mpmdata$TR,
#                  TE = mpmdata$TE,
#                  FA = mpmdata$FA,
#                  TEScale = TEScale,
#                  dataScale = dataScale))
# 
#   ## END function estimateESTATICS()
# }

estimateESTATICSQL <- function (mpmdata, TEScale = 100, dataScale = 1000, sigma = NULL, 
                                L = NULL, verbose = TRUE) 
{
  if (mpmdata$model == 2) {
    xmat <- matrix(0, mpmdata$nFiles, 4)
    xmat[1:length(mpmdata$t1Files), 1] <- 1
    xmat[(length(mpmdata$t1Files) + 1):(length(mpmdata$t1Files) + 
                                          length(mpmdata$mtFiles)), 2] <- 1
    xmat[(length(mpmdata$t1Files) + length(mpmdata$mtFiles) + 
            1):mpmdata$nFiles, 3] <- 1
    xmat[, 4] <- mpmdata$TE/TEScale
    #
    #   check for voxel in mask with all zeros for a modality
    #
    zerovoxel <- as.logical(((apply(mpmdata$ddata[xmat[,1]==1,,,],2:4,sum)==0)|
                               (apply(mpmdata$ddata[xmat[,2]==1,,,],2:4,sum)==0)|
                               (apply(mpmdata$ddata[xmat[,3]==1,,,],2:4,sum)==0))*mpmdata$mask)
  }
  else {
    xmat <- matrix(0, mpmdata$nFiles, 3)
    xmat[1:length(mpmdata$t1Files), 1] <- 1
    xmat[(length(mpmdata$t1Files) + 1):mpmdata$nFiles, 2] <- 1
    xmat[, 3] <- mpmdata$TE/TEScale
    #
    #   check for voxel in mask with all zeros for a modality
    #
    zerovoxel <- as.logical(((apply(mpmdata$ddata[xmat[,1]==1,,,],2:4,sum)==0)|
                               (apply(mpmdata$ddata[xmat[,2]==1,,,],2:4,sum)==0))*mpmdata$mask)
  }
  if (verbose) {
    cat("Design of the model:\n")
    print(xmat)
  }
  
  #
  #   deactivate voxel in mask with all zeros for a modality
  #
  mpmdata$mask[zerovoxel] <- FALSE
  
  sigma <- sigma/dataScale
  CL <- sigma * sqrt(pi/2) * gamma(L + 0.5)/gamma(L)/gamma(1.5)
  R2star <- 0.05 * TEScale
  indT1 <- order(mpmdata$TE[as.logical(xmat[, 1])])[1]
  if (mpmdata$model == 2) {
    indMT <- order(mpmdata$TE[as.logical(xmat[, 2])])[1] + 
      sum(xmat[, 1])
    indPD <- order(mpmdata$TE[as.logical(xmat[, 3])])[1] + 
      sum(xmat[, 1]) + sum(xmat[, 2])
    npar <- 4
  }
  else {
    indPD <- order(mpmdata$TE[as.logical(xmat[, 2])])[1] + 
      sum(xmat[, 1])
    npar <- 3
  }
  if(length(sigma)==1){
      homsigma <- TRUE
      sig <- sigma
  } else if(all(dim(sigma)==mpmdata$sdim)) {
      homsigma <- FALSE
  } else {
      stop("wrong dimension of array sigma")
  }
  isConv <- array(FALSE, mpmdata$sdim)
  modelCoeff <- array(0, c(npar, mpmdata$sdim))
  invCov <- array(0, c(npar, npar, mpmdata$sdim))
  if (verbose) 
    Sys.time()
  for (z in 1:mpmdata$sdim[3]) {
    for (y in 1:mpmdata$sdim[2]) {
      for (x in 1:mpmdata$sdim[1]) {
        if (mpmdata$mask[x, y, z]) {
          if(!homsigma) sig <- sigma[x, y, z]
          ivec <- mpmdata$ddata[, x, y, z]/dataScale
          if (mpmdata$model == 2) {
            th <- c(ivec[indT1] * exp(-xmat[indT1, 4] * 
                                        R2star), ivec[indMT] * exp(-xmat[indMT, 
                                                                         4] * R2star), ivec[indPD] * exp(-xmat[indPD, 
                                                                                                               4] * R2star), R2star)
            res <- try(nls(ivec ~ qflashplQL(par, xmat, 
                                             CL, sig, L), start = list(par = th), 
                           control = list(maxiter = 200, warnOnly = TRUE)))
            if (class(res) == "try-error" || !res$convInfo$isConv || 
                any(coefficients(res) < 0)) 
              #                        cat("nls failed at x",x,"y",y,"z",z,"\n") 
              #                    print(res)
              res <- nls(ivec ~ qflashplQL(par, xmat, 
                                           CL, sig, L), start = list(par = th), 
                         algorithm = "port", control = list(warnOnly = TRUE, 
                                                            printEval = TRUE), lower = rep(0, 4))
            #                    cat("result of port")
            #                    print(res)
          }
          else {
            th <- c(ivec[indT1] * exp(-xmat[indT1, 3] * 
                                        R2star), ivec[indPD] * exp(-xmat[indPD, 
                                                                         3] * R2star), R2star)
            res <- try(nls(ivec ~ qflashpl2QL(par, xmat, 
                                              CL, sig, L), start = list(par = th), 
                           control = list(maxiter = 200, warnOnly = TRUE)))
            if (class(res) == "try-error" || !res$convInfo$isConv || 
                any(coefficients(res) < 0)) 
              #                        cat("nls failed at x",x,"y",y,"z",z,"\n") 
              #                      print(res)
              res <- nls(ivec ~ qflashpl2QL(par, xmat, 
                                            CL, sig, L), start = list(par = th), 
                         algorithm = "port", control = list(warnOnly = TRUE, 
                                                            printEval = TRUE), lower = rep(0, 3))
            #                      cat("result of port")
            #                      print(res)
          }
          sres <- try(summary(res))
          if(class(sres)!="try-error"){
            isConv[x, y, z] <- res$convInfo$isConv
            modelCoeff[, x, y, z] <- sres$coefficients[, 
                                                       1]
            if (sres$sigma != 0) {
              invCovtmp <- try(solve(sres$cov.unscaled))
              if(class(invCovtmp) != "try-error"){
                invCov[, , x, y, z] <- invCovtmp/sres$sigma^2
              } 
            } # else {
            #                      cat("solve(sres$cov.unscaled) failed at x",x,"y",y,"z",z,"\n")
            #                      print(sres$cov.unscaled)
            #                      mpmdata$mask[x,y,z] <- FALSE
            #                  }
          } else {
            mpmdata$mask[x,y,z] <- FALSE
          }
        }
      }
    }
    if (verbose) 
      cat(z, format(Sys.time()), "\n")
  }
  if (verbose) 
    Sys.time()
  invisible(list(modelCoeff = modelCoeff, invCov = invCov, 
                 isConv = isConv, sdim = mpmdata$sdim, nFiles = mpmdata$nFiles, 
                 t1Files = mpmdata$t1Files, pdFiles = mpmdata$pdFiles, 
                 mtFiles = mpmdata$mtFiles, model = mpmdata$model, maskFile = mpmdata$maskFile, 
                 mask = mpmdata$mask, sigma = sigma, L = L, TR = mpmdata$TR, TE = mpmdata$TE, 
                 FA = mpmdata$FA, TEScale = TEScale, dataScale = dataScale))
}

smoothESTATICS <- function(mpmESTATICSModel,
                           mpmData = NULL,
                           kstar = 16,
                           alpha = 0.05,
                           wghts = NULL,
                           verbose = TRUE) {

  ## length of the vector to smooth (# parameters of model)
  nv <- if (mpmESTATICSModel$model == 2) 4 else 3

  ## determine a suitable adaptation bandwidth
  lambda <- nv * qf(1 - alpha, nv, mpmESTATICSModel$nFiles - nv)

  ## adjust for non-isotropic voxel if necessary
  if (is.null(wghts)) wghts <- c(1, 1, 1)
  ## make first spatial dimension unit, and use second and third for reference
  wghts <- wghts[1] / wghts[2:3]

  ## spatial dimension and number of voxel
  n1 <- mpmESTATICSModel$sdim[1]
  n2 <- mpmESTATICSModel$sdim[2]
  n3 <- mpmESTATICSModel$sdim[3]
  n <- n1 * n2 * n3

  ## initialization for first step
  zobj <- list(bi = rep(1, n), theta = mpmESTATICSModel$modelCoeff)
  bi <- zobj$bi

  ## find the number of usable cores
  mc.cores <- setCores(, reprt = FALSE)

  ## define the maximum bandwith from the number of iterations
  hmax <- 1.25^(kstar/3)

  ## for the verbose mode we compute MAD from original data
  if (verbose) {
    mae <- NULL
    protocol <- matrix("", kstar, 1, dimnames = list(paste("step", 1:kstar), "protocol"))
  }

  ## perform the iteration
  k <- 1
  if (verbose) pb <- txtProgressBar(min = 0, max = kstar, style = 3)

  ## if we plan to smooth the original data too, we take special care in the last step
  if (!is.null(mpmData)) {
    if(length(dim(mpmData$ddata)) != 4 | any(mpmESTATICSModel$sdim != dim(mpmData$ddata)[-1]))
      stop("incompatible dimensions of model parameters and original data")
    kstar <- kstar - 1
    smoothData <- TRUE
  } else {
    smoothData <- FALSE
  }

  while (k <= kstar) {
    ## determine the actual bandwidth for this step
    hakt <- gethani(1, 1.25*hmax, 2, 1.25^k, wghts, 1e-4)

    ## we need the (approx.) size of the weigthing scheme array
    dlw <- (2*trunc(hakt/c(1, wghts))+1)[1:3]

    ## perform the actual adaptive smoothing
    zobj <- .Fortran("vaws",
                     as.double(mpmESTATICSModel$modelCoeff),
                     as.logical(mpmESTATICSModel$mask),
                     as.integer(nv),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     hakt = as.double(hakt),
                     as.double(lambda),
                     as.double(zobj$theta),
                     as.double(mpmESTATICSModel$invCov),
                     bi = as.double(zobj$bi),
                     theta = double(nv*n),
                     as.integer(mc.cores),
                     double(prod(dlw)),
                     as.double(wghts),
                     double(nv * mc.cores),
                     PACKAGE = "qMRI")[c("bi", "theta", "hakt")]

    ## use maximum ni
    bi <- zobj$bi <- pmax(bi, zobj$bi)

    ## some verbose stuff
    if (verbose) {
      protocol[k] <- paste("bandwidth: ", signif(hakt, 3),
                           "MSE: ", signif(mean((zobj$theta - mpmESTATICSModel$modelCoeff)^2),3),
                           "MAE: ", m1 <- signif(mean(abs(zobj$theta - mpmESTATICSModel$modelCoeff)),3),
                           "mean(bi):", signif(mean(zobj$bi),3))
      mae <- c(mae, m1)
      setTxtProgressBar(pb, k)
    }

    ## go for next iteration
    k <- k+1
    gc()

  }
  if (smoothData) { ##  modified last step; smoothing data, too

    ## we need the number of files for some array dimensions
    nve <- mpmESTATICSModel$nFiles

    ## determine the actual bandwidth for this step
    hakt <- gethani(1, 1.25*hmax, 2, 1.25^k, wghts, 1e-4)

    ## we need the (approx.) size of the weigthing scheme array
    dlw <- (2*trunc(hakt/c(1, wghts))+1)[1:3]

    ## perform the actual adaptive smoothing
    zobj <- .Fortran("vawsext",
                     as.double(mpmESTATICSModel$modelCoeff),
                     as.logical(mpmESTATICSModel$mask),
                     as.integer(nv),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     as.double(mpmData$ddata),
                     as.integer(nve),
                     hakt = as.double(hakt),
                     as.double(lambda),
                     as.double(zobj$theta),
                     as.double(mpmESTATICSModel$invCov),
                     bi = as.double(zobj$bi),
                     theta = double(nv*n),
                     thext = double(nve*n),
                     as.integer(mc.cores),
                     double(prod(dlw)),
                     as.double(wghts),
                     double(nv*mc.cores),
                     double(nve*mc.cores),
                     PACKAGE = "qMRI")[c("bi", "theta", "thext", "hakt")]

    ## assign the smoothed data
    dim(zobj$thext) <- c(mpmESTATICSModel$nFiles, mpmESTATICSModel$sdim)

    mpmDataSmoothed <- zobj$thext

    ## use maximum ni
    zobj$bi <- pmax(bi, zobj$bi)

    ## some verbose stuff
    if (verbose) {
      protocol[k] <- paste("bandwidth: ", signif(hakt, 3),
                           "MSE: ", signif(mean((zobj$theta - mpmESTATICSModel$modelCoeff)^2),3),
                           "MAE: ", m1 <- signif(mean(abs(zobj$theta - mpmESTATICSModel$modelCoeff)),3),
                           "mean(bi):", signif(mean(zobj$bi),3))
      mae <- c(mae, m1)
      setTxtProgressBar(pb, kstar)
    }

  } else {
    mpmDataSmoothed <- NULL
  }
  if (verbose) close(pb)
  if (verbose) print(protocol)

  ## set correct dimensions
  dim(zobj$theta) <- c(nv, n1, n2, n3)
  dim(zobj$bi) <- c(n1, n2, n3)

  ## assign values
  invisible(list(modelCoeff = zobj$theta,
                 invCov = mpmESTATICSModel$invCov,
                 isConv = mpmESTATICSModel$isConv,
                 bi = zobj$bi,
                 smoothPar = c(lambda, hakt, alpha),
                 smoothedData = mpmDataSmoothed,
                 sdim = mpmESTATICSModel$sdim,
                 nFiles = mpmESTATICSModel$nFiles,
                 t1Files = mpmESTATICSModel$t1Files,
                 pdFiles = mpmESTATICSModel$pdFiles,
                 mtFiles = mpmESTATICSModel$mtFiles,
                 model = mpmESTATICSModel$model,
                 maskFile = mpmESTATICSModel$maskFile,
                 mask = mpmESTATICSModel$mask,
                 TR = mpmESTATICSModel$TR,
                 TE = mpmESTATICSModel$TE,
                 FA = mpmESTATICSModel$FA,
                 TEScale = mpmESTATICSModel$TEScale,
                 dataScale = mpmESTATICSModel$dataScale))
  ## END function smoothESTATICS()
}

calculateQI <- function(mpmESTATICSModel,
                        b1File = NULL,
                        TR2 = 0,
                        verbose = TRUE) {

  ## read B1 correction field
  if (!is.null(b1File)) {
    if (verbose) cat("reading B1 correction file from", b1File, "\n")
    b1Map <- readNIfTI(b1File, reorient = FALSE)/100
    b1Map[b1Map < 0] <- 0
    if (any(dim(b1Map) != mpmESTATICSModel$sdim)) stop("dimension of B1 map does not match data dimension")
  } else {
    if (verbose) cat("no B1 correction\n")
    b1Map <- array(1, mpmESTATICSModel$sdim)
  }

  ## get correct flip angles and TR times
  t1FA <- mpmESTATICSModel$FA[1]
  pdFA <- mpmESTATICSModel$FA[length(mpmESTATICSModel$t1Files) + length(mpmESTATICSModel$mtFiles) + 1]
  t1TR <- mpmESTATICSModel$TR[1]
  pdTR <- mpmESTATICSModel$TR[length(mpmESTATICSModel$t1Files) + length(mpmESTATICSModel$mtFiles) + 1]

  ## calculate E1
  if (verbose) cat("calculating R1 ... ")
  alphat1 <- b1Map * t1FA / 180 * pi
  alphapd <- b1Map * pdFA / 180 * pi
  SINalphat1 <- sin(alphat1)
  COSalphat1 <- cos(alphat1)
  SINalphapd <- sin(alphapd)
  COSalphapd <- cos(alphapd)
  rm(alphat1, alphapd)
  if (mpmESTATICSModel$model == 2) {
    enum <- mpmESTATICSModel$modelCoeff[1, , , ] - SINalphat1/SINalphapd * mpmESTATICSModel$modelCoeff[3, , , ]
    denom <- mpmESTATICSModel$modelCoeff[1, , , ] * COSalphat1 - SINalphat1/SINalphapd * mpmESTATICSModel$modelCoeff[3, , , ] * COSalphapd
  } else {
    enum <- mpmESTATICSModel$modelCoeff[1, , , ] - SINalphat1/SINalphapd * mpmESTATICSModel$modelCoeff[2, , , ]
    denom <- mpmESTATICSModel$modelCoeff[1, , , ] * COSalphat1 - SINalphat1/SINalphapd * mpmESTATICSModel$modelCoeff[2, , , ] * COSalphapd
  }
  E1 <- enum/denom
  rm(enum, denom, COSalphapd, SINalphapd)
  R1 <- -log(E1)/t1TR

  ### RF spoiling correction Preibisch and Deichmann MRM 61 (2009) 125-135
  ### These coefficients depend on the sequence!! See getPolynomsP2_ab and
  ### MTprot in VBQ
  P2_a = getPolynomsP2_ab(pdTR, t1TR, pdFA, t1FA, verbose)$P2_a
  P2_b = getPolynomsP2_ab(pdTR, t1TR, pdFA, t1FA, verbose)$P2_b
  R1 = R1 / ((P2_a[1]*b1Map^2 + P2_a[2]*b1Map + P2_a[3]) * R1 + (P2_b[1]*b1Map^2 + P2_b[2]*b1Map + P2_b[3]))
  E1 = exp(- R1 * t1TR)
  ### END spoiling correction

  if (verbose) cat("done\n")

  ## calculate PD
  if (verbose) cat("calculating PD ... ")
  enum <- (1 - COSalphat1 * E1) * mpmESTATICSModel$modelCoeff[1, , , ] * mpmESTATICSModel$dataScale
  denom <- SINalphat1 * (1 - E1)
  PD <- enum/denom
  rm(enum, denom, SINalphat1)
  if (verbose) cat("done\n")

  ## calculate delta
  if (mpmESTATICSModel$model == 2) {
    if (verbose) cat("calculating MT ... ")
    mtFA <- mpmESTATICSModel$FA[length(mpmESTATICSModel$t1Files) + 1]
    mtTR <- mpmESTATICSModel$TR[length(mpmESTATICSModel$t1Files) + 1]
    alphamt <- b1Map * mtFA / 180 * pi
    E1mt <- E1^(mtTR/t1TR)
    E2mt <- E1^(TR2/t1TR)
    enom <- mpmESTATICSModel$modelCoeff[2, , , ]  * mpmESTATICSModel$dataScale - (1 - E2mt) * sin(alphamt) * PD
    denom <- mpmESTATICSModel$modelCoeff[2, , , ]  * mpmESTATICSModel$dataScale * cos(alphamt) *E1mt + PD * (E2mt  - E1mt) * sin(alphamt)
    delta <- 1 - enom / denom
    rm(alphamt, enom, denom)

    ### correction for MT saturation pulse. see Helms ISMRM 23 (2015) 3360
    delta = 100 * delta * (1 - 0.4) / (1 - 0.4 * b1Map) / b1Map^2;

    if (verbose) cat("done\n")
  } else {
    delta <- NULL
  }

  invisible(list(b1Map = b1Map,
                 R1 = R1 * 1000,
                 R2star = if (mpmESTATICSModel$model == 2) 1000 * mpmESTATICSModel$modelCoeff[4, , , ]/mpmESTATICSModel$TEScale else 1000 * mpmESTATICSModel$modelCoeff[3, , , ]/mpmESTATICSModel$TEScale,
                 PD = PD,
                 delta = delta,
                 model = mpmESTATICSModel$model,
                 t1Files = mpmESTATICSModel$t1Files,
                 mtFiles = mpmESTATICSModel$mtFiles,
                 pdFiles = mpmESTATICSModel$pdFiles,
                 mask = mpmESTATICSModel$mask))
}

imageQI <- function(qi,
                    view = 1,
                    slice = 1) {

  mask <- switch(view,
                 qi$mask[slice, , ],
                 qi$mask[, slice, ],
                 qi$mask[, , slice])
  r2star <- switch(view,
                   qi$R2star[slice, , ],
                   qi$R2star[, slice, ],
                   qi$R2star[, , slice])
  r2star[!mask] <- 0
  r1 <- switch(view,
               qi$R1[slice, , ],
               qi$R1[, slice, ],
               qi$R1[, , slice])
  r1[!mask] <- 0
  pd <- switch(view,
               qi$PD[slice, , ],
               qi$PD[, slice, ],
               qi$PD[, , slice])
  pd[!mask] <- 0
  if (qi$model == 2) {
    delta <- switch(view,
                    qi$delta[slice, , ],
                    qi$delta[, slice, ],
                    qi$delta[, , slice])
    delta[!mask] <- 0
  }
  indx <- 1:dim(r2star)[1]
  indy <- 1:dim(r2star)[2]

  if (qi$model == 2) {
    def.par <- par(mfrow = c(2, 2), mar = c(3, 3, 3, 0))
    rimage(indx, indy, r2star, zlim = c(0, 0.05), main = "R2star")
    rimage(indx, indy, r1, zlim = c(0.0002, 0.0015), main = "R1")
    rimage(indx, indy, pd, zlim = c(0, 10000), main = "PD")
    rimage(indx, indy, delta, zlim = c(0, 0.03), main = "MT")
  } else {
    def.par <- par(mfrow = c(2, 2), mar = c(3, 3, 3, 0))
    rimage(indx, indy, r2star, zlim = c(0, 0.05), main = "R2star")
    rimage(indx, indy, r1, zlim = c(0.0002, 0.0015), main = "R1")
    rimage(indx, indy, pd, zlim = c(0, 10000), main = "PD")
  }



  par(def.par)
}


writeQI <- function(qi,
                    dir = NULL,
                    verbose = TRUE) {

  if (!is.null(dir)) {
    dir.create(dir)
    r2file <- file.path(dir, "R2")
    r1file <- file.path(dir, "R1")
    pdfile <- file.path(dir, "PD")
    mtfile <- file.path(dir, "MT")
  } else {
    r2file <- "R2"
    r1file <- "R1"
    pdfile <- "PD"
    mtfile <- "MT"
  }

  ds <- readNIfTI(qi$t1Files[1], reorient = FALSE)
  ds@datatype <- 16
  ds@magic <- "n+1"
  ds@vox_offset <- 352

  if (verbose) cat("writing R2 file ... ")
  ds@descrip <- "R2"
  writeNIfTI(as.nifti(qi$R2star, ds), file = r2file)
  if (verbose) cat("done\n")
  if (verbose) cat("writing R1 file ... ")
  ds@descrip <- "R1"
  writeNIfTI(as.nifti(qi$R1, ds), file = r1file)
  if (verbose) cat("done\n")
  if (verbose) cat("writing PD file ... ")
  ds@descrip <- "PD"
  writeNIfTI(as.nifti(qi$PD, ds), file = pdfile)
  if (verbose) cat("done\n")
  if (qi$model == 2) {
    if (verbose) cat("writing MT file ... ")
    ds@descrip <- "MT"
    writeNIfTI(as.nifti(qi$delta, ds), file = mtfile)
    if (verbose) cat("done\n")
  }
}

writeESTATICS <- function(mpmESTATICSModel,
                          dir = NULL,
                          prefix = "sm",
                          verbose = TRUE) {

  if (!is.null(dir)) {
    dir.create(dir)
    r2file <- file.path(dir, "R2")
    st1file <- file.path(dir, "ST1")
    spdfile <- file.path(dir, "SPD")
    smtfile <- file.path(dir, "SMT")
  } else {
    r2file <- "R2"
    st1file <- "ST1"
    spdfile <- "SPD"
    smtfile <- "SMT"
  }

  ds <- readNIfTI(mpmESTATICSModel$t1Files[1], reorient = FALSE)
  ds@datatype <- 16
  ds@magic <- "n+1"
  ds@vox_offset <- 352

  if (mpmESTATICSModel$model == 2) {
    if (verbose) cat("writing R2 file ... ")
    ds@descrip <- "R2"
    writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[4, , , ], ds), file = r2file)
    if (verbose) cat("done\n")
    if (verbose) cat("writing ST1 file ... ")
    ds@descrip <- "ST1"
    writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[1, , , ], ds), file = st1file)
    if (verbose) cat("done\n")
    if (verbose) cat("writing SPD file ... ")
    ds@descrip <- "SPD"
    writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[3, , , ], ds), file = spdfile)
    if (verbose) cat("done\n")
    if (verbose) cat("writing SMT file ... ")
    ds@descrip <- "SMT"
    writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[2, , , ], ds), file = smtfile)
    if (verbose) cat("done\n")
  } else {
    if (verbose) cat("writing R2 file ... ")
    ds@descrip <- "R2"
    writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[3, , , ], ds), file = r2file)
    if (verbose) cat("done\n")
    if (verbose) cat("writing ST1 file ... ")
    ds@descrip <- "ST1"
    writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[1, , , ], ds), file = st1file)
    if (verbose) cat("done\n")
    if (verbose) cat("writing SPD file ... ")
    ds@descrip <- "SPD"
    writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[2, , , ], ds), file = spdfile)
    if (verbose) cat("done\n")
  }

  if (!is.null(mpmESTATICSModel$smoothedData)) {
    ii <- 1
    for (i in 1:length(t1Files)) {
      ds <- readNIfTI(t1Files[i], reorient = FALSE)
      ds@datatype <- 16
      ds@magic <- "n+1"
      ds@vox_offset <- 352
      if (!is.null(dir)) {
        fname <- file.path(dir, paste(prefix, basename(t1Files[i]), sep = ""))
      } else {
        fname <- paste(prefix, basename(t1Files[i]), sep = "")
      }
      if (verbose) cat("writing", fname, "... ")
      writeNIfTI(as.nifti(mpmESTATICSModel$smoothedData[ii, , , ], ds), file = fname)
      if (verbose) cat("done\n")
      ii <- ii + 1
    }
    if (mpmESTATICSModel$model == 2) {
      for (i in 1:length(mtFiles)) {
        ds <- readNIfTI(mtFiles[i], reorient = FALSE)
        ds@datatype <- 16
        ds@magic <- "n+1"
        ds@vox_offset <- 352
        if (!is.null(dir)) {
          fname <- file.path(dir, paste(prefix, basename(mtFiles[i]), sep = ""))
        } else {
          fname <- paste(prefix, basename(mtFiles[i]), sep = "")
        }
        if (verbose) cat("writing", fname, "... ")
        writeNIfTI(as.nifti(mpmESTATICSModel$smoothedData[ii, , , ], ds), file = fname)
        if (verbose) cat("done\n")
        ii <- ii + 1
      }
    }
    for (i in 1:length(pdFiles)) {
      ds <- readNIfTI(pdFiles[i], reorient = FALSE)
      ds@datatype <- 16
      ds@magic <- "n+1"
      ds@vox_offset <- 352
      if (!is.null(dir)) {
        fname <- file.path(dir, paste(prefix, basename(pdFiles[i]), sep = ""))
      } else {
        fname <- paste(prefix, basename(pdFiles[i]), sep = "")
      }
      if (verbose) cat("writing", fname, "... ")
      writeNIfTI(as.nifti(mpmESTATICSModel$smoothedData[ii, , , ], ds), file = fname)
      if (verbose) cat("done\n")
      ii <- ii + 1
    }
  }

}

estimateQIconf <- function(mpmESTATICSmodel,
                           verbose = TRUE) {

  if (mpmESTATICSmodel$model != 2) stop("only full model implemented!")

  ## TODO: this should be done in estimateESTATICS and smoothESTATICS
  dimnames(mpmESTATICSmodel$modelCoeff) <- list(c("ST1", "SMT", "SPD", "R2star"), NULL, NULL, NULL)
  dimnames(mpmESTATICSmodel$invCov) <- list(c("ST1", "SMT", "SPD", "R2star"), c("ST1", "SMT", "SPD", "R2star"), NULL, NULL, NULL)

  R1 <- array(0, mpmESTATICSmodel$sdim)
  CIR1 <- array(0, c(2, mpmESTATICSmodel$sdim))
  R2 <- array(0, mpmESTATICSmodel$sdim)
  CIR2 <- array(0, c(2, mpmESTATICSmodel$sdim))
  for (z in 1:mpmESTATICSmodel$sdim[3]) {
    for (y in 1:mpmESTATICSmodel$sdim[2]) {
      for (x in 1:mpmESTATICSmodel$sdim[1]) {
        if (mpmESTATICSmodel$mask[x, y, z]) {
          if (is.null(mpmESTATICSmodel$bi)) {
            zz <- ESTATICS.confidence(mpmESTATICSmodel$modelCoeff[, x, y, z], mpmESTATICSmodel$invCov[, , x, y, z], mpmESTATICSmodel$FA[1]*pi/180, mpmESTATICSmodel$FA[length(mpmESTATICSmodel$t1Files) + length(mpmESTATICSmodel$mtFiles) + 1]*pi/180, mpmESTATICSmodel$TR[1], df=mpmESTATICSmodel$nFiles-4, 0.05)
          } else {
            zz <- ESTATICS.confidence(mpmESTATICSmodel$modelCoeff[, x, y, z], mpmESTATICSmodel$invCov[, , x, y, z] * mpmESTATICSmodel$bi[x, y, z], mpmESTATICSmodel$FA[1]*pi/180, mpmESTATICSmodel$FA[length(mpmESTATICSmodel$t1Files) + length(mpmESTATICSmodel$mtFiles) + 1]*pi/180, mpmESTATICSmodel$TR[1], df=mpmESTATICSmodel$nFiles-4, 0.05)
          }
          R2[x, y, z] <- zz$R2star
          R1[x, y, z] <- zz$R1
          CIR1[, x, y, z] <- zz$CIR1
          CIR2[, x, y, z] <- zz$CIR2star
        }
      }
    }
    if (verbose) cat(z, format(Sys.time()), "\n")
  }

  invisible(list(R1 = R1,
                 CIR1 = CIR1,
                 R2star = R2,
                 CIR2star = CIR2,
                 model = mpmESTATICSmodel$model,
                 t1Files = mpmESTATICSmodel$t1Files,
                 mtFiles = mpmESTATICSmodel$mtFiles,
                 pdFiles = mpmESTATICSmodel$pdFiles,
                 mask = mpmESTATICSmodel$mask))
}

writeQIconf <- function(qiConf,
                        dir= NULL,
                        verbose = TRUE) {
  if (!is.null(dir)) {
    dir.create(dir)
    r1file <- file.path(dir, "R1")
    r1Lfile <- file.path(dir, "R1lower")
    r1Ufile <- file.path(dir, "R1upper")
    r2file <- file.path(dir, "R2")
    r2Lfile <- file.path(dir, "R2lower")
    r2Ufile <- file.path(dir, "R2upper")
  } else {
    r1file <- "R1"
    r1Lfile <- "R1lower"
    r1Ufile <- "R1upper"
    r2file <- "R2"
    r2Lfile <- "R2lower"
    r2Ufile <- "R2upper"
  }

  ds <- readNIfTI(qiConf$t1Files[1], reorient = FALSE)
  ds@datatype <- 16
  ds@magic <- "n+1"
  ds@vox_offset <- 352

  if (verbose) cat("writing R1 file ... ")
  ds@descrip <- "R1"
  writeNIfTI(as.nifti(qiConf$R1, ds), file = r1file)
  if (verbose) cat("done\n")
  if (verbose) cat("writing R1lower file ... ")
  ds@descrip <- "R1lower"
  writeNIfTI(as.nifti(qiConf$CIR1[1, , , ], ds), file = r1Lfile)
  if (verbose) cat("done\n")
  if (verbose) cat("writing R1upper file ... ")
  ds@descrip <- "R1upper"
  writeNIfTI(as.nifti(qiConf$CIR1[2, , , ], ds), file = r1Ufile)
  if (verbose) cat("done\n")

  if (verbose) cat("writing R2 file ... ")
  ds@descrip <- "R2"
  writeNIfTI(as.nifti(qiConf$R2star, ds), file = r2file)
  if (verbose) cat("done\n")
  if (verbose) cat("writing R2lower file ... ")
  ds@descrip <- "R2lower"
  writeNIfTI(as.nifti(qiConf$CIR2star[1, , , ], ds), file = r2Lfile)
  if (verbose) cat("done\n")
  if (verbose) cat("writing R2upper file ... ")
  ds@descrip <- "R2upper"
  writeNIfTI(as.nifti(qiConf$CIR2star[2, , , ], ds), file = r2Ufile)
  if (verbose) cat("done\n")

}

getPolynomsP2_ab <- function(TR_pdw, TR_t1w, fa_pdw, fa_t1w, verbose = TRUE) {

  ## Settings for R. Deichmann steady state correction using T2 = 64ms at 3T
  ## Correction parameters were calculated for 3 different parameter sets:
  if ((TR_pdw == 23.7) && (TR_t1w == 18.7) && (fa_pdw == 6) && (fa_t1w == 20)) {
    ## 1) classic FIL protocol (Weiskopf et al., Neuroimage 2011):
    ## PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=18.7ms; a=20deg
    if (verbose) cat("Spoiling correction ... Classic FIL protocol\n")
    P2_a <- c(78.9228195006542,   -101.113338489192,    47.8783287525126)
    P2_b <- c(-0.147476233142129,    0.126487385091045,  0.956824374979504)
  } else if ((TR_pdw == 24.5) && (TR_t1w == 24.5) && (fa_pdw == 5) && (fa_t1w == 29)) {
    ## 2) new FIL/Helms protocol
    ## PD-weighted: TR=24.5ms; a=5deg; T1-weighted: TR=24.5ms; a=29deg
    if (verbose) cat("Spoiling correction ... New FIL/Helms protocol\n")
    P2_a <- c(93.455034845930480, -120.5752858196904,   55.911077913369060)
    P2_b <- c(-0.167301931434861,    0.113507432776106,  0.961765216743606)
  } else if ((TR_pdw == 24.0) && (TR_t1w == 19.0) && (fa_pdw == 6) && (fa_t1w == 20)) {
    ## 3) Siemens product sequence protocol used in Lausanne (G Krueger)
    ## PD-weighted: TR=24ms; a=6deg; T1-weighted: TR=19ms; a=20deg
    if (verbose) cat("Spoiling correction ... Siemens product Lausanne (GK) protocol\n")
    P2_a <- c(67.023102027100880, -86.834117103841540, 43.815818592349870)
    P2_b <- c(-0.130876849571103,   0.117721807209409,  0.959180058389875)
  } else if ((TR_pdw == 23.7) && (TR_t1w == 23.7) && (fa_pdw == 6) && (fa_t1w == 28)) {
    ## 4) High-res (0.8mm) FIL protocol:
    ## PD-weighted: TR=23.7ms; a=6deg; T1-weighted: TR=23.7ms; a=28deg
    if (verbose) cat("Spoiling correction ... High-res FIL protocol\n")
    P2_a <- c( 1.317257319014170e+02, -1.699833074433892e+02, 73.372595677371650)
    P2_b <- c(-0.218804328507184,      0.178745853134922,      0.939514554747592)
  } else if ((TR_pdw == 25.25) && (TR_t1w == 25.25) && (fa_pdw == 5) && (fa_t1w == 29)) {
    ## 4)NEW  High-res (0.8mm) FIL protocol:
    ## PD-weighted: TR=25.25ms; a=5deg; T1-weighted: TR=TR=25.25ms; a=29deg
    if (verbose) cat("Spoiling correction ... High-res FIL protocol\n")
    P2_a <- c(88.8623036106612,   -114.526218941363,    53.8168602253166)
    P2_b <- c(-0.132904017579521,    0.113959390779008,  0.960799295622202)
  } else if ((TR_pdw == 24.5) && (TR_t1w == 24.5) && (fa_pdw == 6) && (fa_t1w == 21)) {
    ## 5)NEW  1mm protocol - seq version v2k:
    ## PD-weighted: TR=24.5ms; a=6deg; T1-weighted: TR=24.5ms; a=21deg
    if (verbose) cat("Spoiling correction ... v2k protocol")
    P2_a <- c(71.2817617982844,   -92.2992876164017,   45.8278193851731)
    P2_b <- c(-0.137859046784839,   0.122423212397157,  0.957642744668469)
  } else if ((TR_pdw == 25.0) && (TR_t1w == 25.0) && (fa_pdw == 6) && (fa_t1w == 21)) {
    ## 6) 800um protocol - seq version v3* released used by MEG group:
    ## TR = 25ms for all volumes; flipAngles = [6, 21 deg] for PDw and T1w
    ## Correction parameters below were determined via Bloch-Torrey
    ## simulations but end result agrees well with EPG-derived correction
    ## for this RF spoiling increment of 137 degrees.
    ## See: Callaghan et al. ISMRM, 2015, #1694
    if (verbose) cat("Spoiling correction ... v3* 0.8mm R4 protocol\n")
    P2_a <- c(57.427573706259864, -79.300742898810441,  39.218584751863879)
    P2_b <- c(-0.121114060111119,   0.121684347499374,   0.955987357483519)
  } else {
    if (verbose) cat("Spoiling correction ... not defined for this protocol. No correction being applied.\n")
    P2_a <- c(0, 0, 0)
    P2_b <- c(0, 0, 1)
  }
  list(P2_a = P2_a, P2_b = P2_b)
}
