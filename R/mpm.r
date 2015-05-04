mpmData <- function(myt1Files  = NULL,
                    mypdFiles  = NULL,
                    mymtFiles  = NULL,
                    mymaskFile = NULL,
                    myb1File   = NULL,
                    mysdim     = NULL,
                    verbose    = TRUE) {
  
  
  ## Get the environment for this
  ## instance of the function.
  thisEnv <- environment()
  
  ## we need at least T1w and PDw files
  if (is.null(myt1Files)) stop("vector of T1 files required")
  if (is.null(mypdFiles)) stop("vector of PD files required")
  ## TODO: test whether there are enough files for the model?

  if (is.null(mysdim)) stop("need spatial dimensionality of the data")
  if (!is.numeric(mysdim) | length(mysdim) != 3) stop("need exactly three numbers for spatial dimensions")
  
  t1Files <- myt1Files
  pdFiles <- mypdFiles
  mtFiles <- mymtFiles
  maskFile <- mymaskFile
  b1File <- myb1File

  ## select the model according to the existence of MTw files
  model <- if (is.null(mymtFiles)) {
    1 # the simple model without MT inclusion
  } else {
    2 # the model including MT
  }
  
  ## count the number of data volumes (if (is.null(mtFiles)) length(mtFiles) == 0)
  nFiles <- length(t1Files) + length(mtFiles) + length(pdFiles)
  
  # this is the spatial dimensionality of the data
  sdim <- mysdim
  
  # the array for the data itself
  ddata <- numeric(nFiles * prod(sdim))
  dim(ddata) <- c(nFiles, prod(sdim))
  ddataSmoothed <- NULL
  
  ## for each files we have a TR, TE, and flip angle (FA)
  TR <- TE <- FA <- numeric(nFiles)
  #   t1TR <- t1TE <- t1FA <- numeric(length(get("t1Files", thisEnv)))
  #   pdTR <- pdTE <- pdFA <- numeric(length(get("pdFiles", thisEnv)))
  #   if (get("model", thisEnv) == 2) {
  #     mtTR <- mtTE <- mtFA <- numeric(length(get("mtFiles", thisEnv)))
  #   }
  
  ## ... now we read all data volumes and extract the TR, TE, and FA values for each ...
  ii <- 1
  ## ... for all T1 volumes ...
  if (verbose) cat("reading T1 files\n")
  if (verbose) pb <- txtProgressBar(min = 0, max = length(t1Files), style = 3)
  for (i in 1:length(t1Files)) {
    ds <- readNIfTI(t1Files[i], reorient = FALSE)
    ddata[ii, ] <- ds
    ## IMPORTANT: This is special to Siawoosh data
    res <- str_match_all(ds@descrip, "([[:alpha:]]{2})=([.0123456789]+)([[:alpha:]]{2,})")
    for (nn in 1:dim(res[[1]])[1]) {
      if (res[[1]][nn, 2] == "TR") TR[ii] <- as.numeric(res[[1]][nn, 3])
      if (res[[1]][nn, 2] == "TE") TE[ii] <- as.numeric(res[[1]][nn, 3])
      if (res[[1]][nn, 2] == "FA") FA[ii] <- as.numeric(res[[1]][nn, 3])
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
      ## IMPORTANT: This is special to Siawoosh data
      res <- str_match_all(ds@descrip, "([[:alpha:]]{2})=([.0123456789]+)([[:alpha:]]{2,})")
      for (nn in 1:dim(res[[1]])[1]) {
        if (res[[1]][nn, 2] == "TR") TR[ii] <- as.numeric(res[[1]][nn, 3])
        if (res[[1]][nn, 2] == "TE") TE[ii] <- as.numeric(res[[1]][nn, 3])
        if (res[[1]][nn, 2] == "FA") FA[ii] <- as.numeric(res[[1]][nn, 3])
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
    ## IMPORTANT: This is special to Siawoosh data
    res <- str_match_all(ds@descrip, "([[:alpha:]]{2})=([.0123456789]+)([[:alpha:]]{2,})")
    for (nn in 1:dim(res[[1]])[1]) {
      if (res[[1]][nn, 2] == "TR") TR[ii] <- as.numeric(res[[1]][nn, 3])
      if (res[[1]][nn, 2] == "TE") TE[ii] <- as.numeric(res[[1]][nn, 3])
      if (res[[1]][nn, 2] == "FA") FA[ii] <- as.numeric(res[[1]][nn, 3])
    }
    ii <- ii + 1
    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) close(pb)
  dim(ddata) <- c(nFiles, sdim)
  ## ... done!  
  
  # the array for the model parameters and inverse covariance
  if (model == 1) {
    modelCoeff <- numeric(3 * prod(sdim))
    dim(modelCoeff) <- c(3, sdim)
    invCov <- numeric(3 * 3 * prod(sdim))
    dim(invCov) <- c(3, 3, sdim)
    pNames <- c("ST1", "SPD", "R2star")
    dimnames(modelCoeff) <- list(pNames, NULL, NULL, NULL)
    dimnames(invCov) <- list(pNames, pNames, NULL, NULL, NULL)
    QI <- numeric(3 * prod(sdim))
    dim(QI) <- c(3, sdim)
    QIsmoothed <- numeric(3 * prod(sdim))
    dim(QIsmoothed) <- c(3, sdim)
  } else {
    modelCoeff <- numeric(4 * prod(sdim))    
    dim(modelCoeff) <- c(4, sdim)
    invCov <- numeric(4 * 4 * prod(sdim))
    dim(invCov) <- c(4, 4, sdim)
    pNames <- c("ST1", "SMT", "SPD", "R2star")
    dimnames(modelCoeff) <- list(pNames, NULL, NULL, NULL)
    dimnames(invCov) <- list(pNames, pNames, NULL, NULL, NULL)
    QI <- numeric(4 * prod(sdim))
    dim(QI) <- c(4, sdim)
    QIsmoothed <- numeric(4 * prod(sdim))
    dim(QIsmoothed) <- c(4, sdim)
  }
  modelCoeffSmoothed <- NULL
  
  ## we need some scales to improve the condition number of the covariance matrix:
  TEScale <- 100
  dataScale <- 1000

  # an array, whether the model estimation converged
  isConv <- array(FALSE, sdim)
  
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
  
  # the b1 File
  ## TODO: set b1Map to 1 also if b1File does not exist or cannot be read by readNIFTI()
  if (is.null(b1File)) {
    if (verbose) cat("found no b1 file, setting mask TRUE everywhere")
    b1Map <- array(1, sdim)
    ## TODO: do we need a full array here, or could a single value be sufficient?
  } else {
    if (verbose) cat("reading B1 file ... ")
    b1Map <- readNIfTI(b1File, reorient = FALSE)/100
    b1Map[b1Map < 0] <- 0    
    if (verbose) cat("done\n")
  }
  
  ## some names for smoothing parameters and results
  lambda <- NULL
  hakt <- NULL
  mae <- NULL
  bi <- NULL
  alpha <- NULL
  
  me <- list(
    
    ## Define the environment where this list is defined so
    ## that I can refer to it later.
    thisEnv = thisEnv,
    
    ## Define the accessors for the data fields.
    getEnv = function()
    {
      return(get("thisEnv", thisEnv))
    },
    
    getT1Files = function()
    {
      return(get("t1Files", thisEnv))
    },
    
    getPDFiles = function()
    {
      return(get("pdFiles", thisEnv))
    },
    
    getMTFiles = function()
    {
      return(get("mtFiles", thisEnv))
    },
    
    getMaskFile = function()
    {
      return(get("maskFile", thisEnv))
    },
    
    setMaskFile = function(mymaskFile) {
      ## TODO: I am not sure, whether we need the assign() function here.
      assign("mask", as.logical(readNIfTI(mymaskFile, reorient = FALSE)), thisEnv)
      ## TODO: I am not sure, whether we need the get() function here.
      dim(get("mask", thisEnv)) <- get("sdim", thisEnv)
      assign("maskFile", mymaskFile, thisEnv)
    },
    
    setMask = function(mymask) {
      assign("mask", mymask, thisEnv)
      dim(get("mask", thisEnv)) <- get("sdim", thisEnv)
      assign("maskFile", NULL, thisEnv)
    },
    
    getB1File = function()
    {
      return(get("b1File", thisEnv))
    },
    
    setB1File = function(myb1File) {
      ## TODO: I am not sure, whether we need the assign() function here.
      assign("b1Map", readNIfTI(myb1File, reorient = FALSE)/100, thisEnv)
      get("b1Map", thisEnv)[get("b1Map", thisEnv) < 0] <- 0
      ## TODO: I am not sure, whether we need the get() function here.
      dim(get("b1Map", thisEnv)) <- get("sdim", thisEnv)
      assign("b1File", myb1File, thisEnv)
    },
    
    setB1Map = function(myb1Map) {
      assign("b1Map", myb1Map, thisEnv)
      dim(get("b1Map", thisEnv)) <- get("sdim", thisEnv)
      assign("b1File", NULL, thisEnv)
    },
    
    getMPMData = function() {
      return(get("ddata", thisEnv))
    },
    
    getTR = function() {
      return(get("TR", thisEnv))
    },
    
    getTE = function() {
      return(get("TE", thisEnv))
    },
    
    getFA = function() {
      return(get("FA", thisEnv))
    },
    
    getSdim = function() {
      return(get("sdim", thisEnv))
    },
    
    getDdim = function() {
      return(c(get("nFiles", thisEnv), get("sdim", thisEnv)))
    },
    
    getScales = function() {
      cat("Scale for TE times   :", get("TEScale", thisEnv), "\n")
      cat("Scale for Data vector:", get("dataScale", thisEnv), "\n")
      return(c(get("TEScale", thisEnv), get("dataScale", thisEnv)))
    },
    
    setScales = function(myTEScale = NULL, mydataScale = NULL) {
      if (!is.null(myTEScale)) assign("TEScale", myTEScale, thisEnv)
      if (!is.null(mydataScale)) assign("dataScale", mydataScale, thisEnv)
    },
    
    estimateESTATICS = function(verbose = TRUE) {
      
      ## this is our design ...
      if (model == 2) {
        xmat <- matrix(0, get("nFiles", thisEnv), 4)
        xmat[1:length(get("t1Files", thisEnv)), 1] <- 1
        xmat[(length(get("t1Files", thisEnv))+1):(length(get("t1Files", thisEnv))+length(get("mtFiles", thisEnv))), 2] <- 1
        xmat[(length(get("t1Files", thisEnv))+length(get("mtFiles", thisEnv))+1):(length(get("t1Files", thisEnv))+length(get("mtFiles", thisEnv))+length(get("pdFiles", thisEnv))), 3] <- 1
        xmat[, 4] <- get("TE", thisEnv) / get("TEScale", thisEnv)
        ## ... for our model in qflashpl() ...
        ## S_{T1} = par[1] * exp(- par[4] * TE)
        ## S_{MT} = par[2] * exp(- par[4] * TE)
        ## S_{PD} = par[3] * exp(- par[4] * TE)
      } else {
        xmat <- matrix(0, get("nFiles", thisEnv), 3)
        xmat[1:length(get("t1Files", thisEnv)), 1] <- 1
        xmat[(length(get("t1Files", thisEnv))+1):(length(get("t1Files", thisEnv))+length(get("pdFiles", thisEnv))), 2] <- 1
        xmat[, 3] <- get("TE", thisEnv) / get("TEScale", thisEnv)
        ## ... for our model in qflashpl2() ...
        ## S_{T1} = par[1] * exp(- par[3] * TE)
        ## S_{PD} = par[2] * exp(- par[3] * TE)        
      }
      if (verbose) {
        cat("Design of the model:\n")
        print(xmat)
      }
        
      ## starting value for R* estimate
      R2star <- 0.05 * get("TEScale", thisEnv)
      indT1 <- order(get("TE", thisEnv)[as.logical(xmat[, 1])])[1]
      if (get("model", thisEnv) == 2) {
        indMT <- order(get("TE", thisEnv)[as.logical(xmat[, 2])])[1] + sum(xmat[, 1])
        indPD <- order(get("TE", thisEnv)[as.logical(xmat[, 3])])[1] + sum(xmat[, 1]) + sum(xmat[, 2]) 
      } else {
        indPD <- order(get("TE", thisEnv)[as.logical(xmat[, 2])])[1] + sum(xmat[, 1])        
      }
      
      ## now perform the voxelwise regression
      if (verbose) Sys.time()
      for (z in 1:get("sdim", thisEnv)[3]){
        for (y in 1:get("sdim", thisEnv)[2]) {
          for (x in 1:get("sdim", thisEnv)[1]) {
            if (mask[x, y, z]) {
              ivec  <- get("ddata", thisEnv)[, x, y, z] / get("dataScale", thisEnv)
              if (get("model", thisEnv) == 2) { 
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
              get("isConv", thisEnv)[x, y, z] <- res$convInfo$isConv
              get("modelCoeff", thisEnv)[, x, y, z] <- sres$coefficients[, 1]
              if (sres$sigma != 0) {
                get("invCov", thisEnv)[, , x, y, z] <- solve(sres$cov.unscaled) / sres$sigma^2
              } else {
                ## TODO: CHECK!
                get("invCov", thisEnv)[, , x, y, z] <- 0
              }
            }
          }
        }
        if (verbose) cat(z, format(Sys.time()), "\n")
      }
      if (verbose) Sys.time()

      ## END function estimateESTATICS()
    },
    
    smoothESTATICS = function(kstar = 16, 
                              smoothData = FALSE,
                              alpha = 0.05, 
                              wghts = NULL,
                              verbose = TRUE) {
      
      ## length of the vector to smooth (# parameters of model)
      nv <- if (get("model", thisEnv) == 2) 4 else 3

      ## determine a suitable adaptation bandwidth
      lambda <- nv * qf(1 - alpha, nv, get("nFiles", thisEnv) - nv)
      
      ## adjust for non-isotropic voxel if necessary
      if (is.null(wghts)) wghts <- c(1, 1, 1)
      ## make first spatial dimension unit, and use second and third for reference
      wghts <- wghts[1] / wghts[2:3]
      
      ## spatial dimension and number of voxel
      n1 <- get("sdim", thisEnv)[1]
      n2 <- get("sdim", thisEnv)[2]
      n3 <- get("sdim", thisEnv)[3]
      n <- n1 * n2 * n3
      
      ## initialization for first step
      zobj <- list(bi = rep(1, n), theta = get("modelCoeff", thisEnv))
      bi <- zobj$bi
      
      ## if we plan to smooth the original data too, we take special care in the last step
      if (smoothData) kstar <- kstar-1
      
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
      while (k <= kstar) {
        ## determine the actual bandwidth for this step
        hakt <- gethani(1, 1.25*hmax, 2, 1.25^k, wghts, 1e-4)
        
        ## we need the (approx.) size of the weigthing scheme array 
        dlw <- (2*trunc(hakt/c(1, wghts))+1)[1:3]

        ## perform the actual adaptive smoothing 
        zobj <- .Fortran("vaws",
                         as.double(get("modelCoeff", thisEnv)),
                         as.logical(get("mask", thisEnv)),
                         as.integer(nv),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(n3),
                         hakt = as.double(hakt),
                         as.double(lambda),
                         as.double(zobj$theta),
                         as.double(get("invCov", thisEnv)),
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
                               "MSE: ", signif(mean((zobj$theta - get("modelCoeff", thisEnv))^2),3),
                               "MAE: ", m1 <- signif(mean(abs(zobj$theta - get("modelCoeff", thisEnv))),3),
                               "mean(bi):", signif(mean(zobj$bi),3),
                               "\n")
          mae <- c(mae, m1)
          setTxtProgressBar(pb, i)
        }

        ## go for next iteration
        k <- k+1
        gc()
        
      }
      if (smoothData) { ##  modified last step; smoothing data, too

        ## we need the number of files for some array dimensions
        nve <- get("nFiles", thisEnv)

        ## determine the actual bandwidth for this step
        hakt <- gethani(1, 1.25*hmax, 2, 1.25^k, wghts, 1e-4)
        
        ## we need the (approx.) size of the weigthing scheme array 
        dlw <- (2*trunc(hakt/c(1, wghts))+1)[1:3]

        ## perform the actual adaptive smoothing
        zobj <- .Fortran("vawsext",
                         as.double(get("modelCoeff", thisEnv)),
                         as.logical(get("mask", thisEnv)),
                         as.integer(nv),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(n3),
                         as.double(get("ddata", thisEnv)),
                         as.integer(nve),
                         hakt = as.double(hakt),
                         as.double(lambda),
                         as.double(zobj$theta),
                         as.double(get("invCov", thisEnv)),
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
        dim(zobj$thext) <- c(get("nFiles", thisEnv), get("sdim", thisEnv))
        assign("ddataSmoothed", zobj$thext, thisEnv)

        ## use maximum ni
        zobj$bi <- pmax(bi, zobj$bi)

        ## some verbose stuff
        if (verbose) {
          protocol[k] <- paste("bandwidth: ", signif(hakt, 3),
                               "MSE: ", signif(mean((zobj$theta - get("modelCoeff", thisEnv))^2),3),
                               "MAE: ", m1 <- signif(mean(abs(zobj$theta - get("modelCoeff", thisEnv))),3),
                               "mean(bi):", signif(mean(zobj$bi),3),
                               "\n")
          mae <- c(mae, m1)
          setTxtProgressBar(pb, i)
        }
        
      } else {
        thest <- NULL
      }
      if (verbose) close(pb)
      if (verbose) print(protocol)
      
      ## set correct dimensions
      dim(zobj$theta) <- c(nv, n1, n2, n3)
      dim(zobj$bi) <- c(n1, n2, n3)

      ## assign values
      assign("lambda", lambda, thisEnv)
      assign("hakt", hakt, thisEnv)
      assign("alpha", alpha, thisEnv)
      assign("mae", mae, thisEnv)
      assign("modelCoeffSmoothed", zobj$theta, thisEnv)      
      assign("bi", zobj$bi, thisEnv)      

      ## END function smoothESTATICS()      
    },
    
    getQMapsFromESTATICS <- function(t1FA = NULL, 
                                     mtFA = NULL,
                                     pdFA = NULL,
                                     t1TR = NULL,
                                     mtTR1 = NULL, 
                                     mtTR2 = NULL) {
      
      cc <- if (get("model", thisEnv) == 2) 3 else 2

      ## calculate f1
      alphat1 <- get("b1Map", thisEnv) * t1FA / 180 * pi
      alphapd <- get("b1Map", thisEnv)  * pdFA / 180 * pi
      SINalphat1 <- sin(alphat1)
      COSalphat1 <- cos(alphat1)
      SINalphapd <- sin(alphapd)
      COSalphapd <- cos(alphapd)
      if (!is.null(lambda)) {
        enum.s <- get("modelCoeffSmoothed", thisEnv)[1, , , ] - SINalphat1/SINalphapd * get("modelCoeffSmoothed", thisEnv)[cc, , , ]
        denom.s <- get("modelCoeffSmoothed", thisEnv)[1, , , ] * COSalphat1 - SINalphat1/SINalphapd * get("modelCoeffSmoothed", thisEnv)[cc, , , ] * COSalphapd
        f1.s <- enum.s/denom.s
        R1.s <- -log(f1.s) / t1TR
      }
      enum <- get("modelCoeff", thisEnv)[1, , , ] - SINalphat1/SINalphapd * get("modelCoeff", thisEnv)[cc, , , ]
      denom <- get("modelCoeff", thisEnv)[1, , , ] * COSalphat1 - SINalphat1/SINalphapd * get("modelCoeff", thisEnv)[cc, , , ] * COSalphapd
      f1 <- enum/denom
      R1 <- -log(f1) / t1TR
      
      ## calculate PD
      if (!is.null(lambda)) {
        enum.s <- (1 - COSalphat1 * f1.s) * get("modelCoeffSmoothed", thisEnv)[1, , , ]
        denom.s <- SINalphat1 * (1 - f1.s)
        PD.s <- enum.s/denom.s
      }
      enum <- (1 - COSalphat1 * f1) * get("modelCoeff", thisEnv)[1, , , ]
      denom <- SINalphat1 * (1 - f1)
      PD <- enum/denom
      
      if (get("model", thisEnv) == 2) {
        ## calculate delta
        alphamt <- get("b1Map", thisEnv)  * mtFA /180 * pi
        if (!is.null(lambda)) {
          f2.s <- f1.s^(mtTR2/mtTR1)
          enom.s <- get("modelCoeffSmoothed", thisEnv)[2, , , ] - (1 - f2.s) * sin(alphamt) * PD.s
          denom.s <- get("modelCoeffSmoothed", thisEnv)[2, , , ] * cos(alphamt) *f1.s*f2.s - PD.s * f2.s * (1 - f1.s) * sin(alphamt)
          delta.s <- enom.s / denom.s
        }
        f2 <- f1^(mtTR2/mtTR1)
        enom <- get("modelCoeffSmoothed", thisEnv)[2, , , ] - (1 - f2) * sin(alphamt) * PD
        denom <- get("modelCoeffSmoothed", thisEnv)[2, , , ] * cos(alphamt) *f1*f2 - PD * f2 * (1 - f1) * sin(alphamt)
        delta <- enom / denom
      }
      
      if (!is.null(lambda)) {
        get("QISmoothed", thisEnv)[1, , , ] <- R1.s 
        get("QISmoothed", thisEnv)[2, , , ] <- PD.s
        if (get("model", thisEnv) == 2) {
          get("QISmoothed", thisEnv)[3, , , ] <- get("modelCoeffSmoothed", thisEnv)[4, , , ]
          get("QISmoothed", thisEnv)[4, , , ] <- delta.s
        } else {
          get("QISmoothed", thisEnv)[3, , , ] <- get("modelCoeffSmoothed", thisEnv)[3, , , ]          
        }
      } else {
        get("QI", thisEnv)[1, , , ] <- R1.s 
        get("QI", thisEnv)[2, , , ] <- PD.s
        if (get("model", thisEnv) == 2) {
          get("QI", thisEnv)[3, , , ] <- get("modelCoeff", thisEnv)[4, , , ]
          get("QI", thisEnv)[4, , , ] <- delta.s
        } else {
          get("QI", thisEnv)[3, , , ] <- get("modelCoeff", thisEnv)[3, , , ]          
        }        
      }
          
    },
    
    writeQMaps = function() {
      for (i in 1:length(t1Files)) {
        ds <- readNIfTI(t1Files[i], reorient = FALSE)
        ds@datatype <- 16
        writeNIfTI(as.nifti(rcoef.s$yextsmoothed[i, , , ], ds), file = paste(t1Files[i], "sm", sep=""))
      }
      for (i in 1:length(mtFiles)) {
        ds <- readNIfTI(mtFiles[i], reorient = FALSE)
        ds@datatype <- 16
        writeNIfTI(as.nifti(rcoef.s$yextsmoothed[i+8, , , ], ds), file = paste(mtFiles[i], "sm", sep=""))
      }
      for (i in 1:length(pdFiles)) {
        ds <- readNIfTI(pdFiles[i], reorient = FALSE)
        ds@datatype <- 16
        writeNIfTI(as.nifti(rcoef.s$yextsmoothed[i+14, , , ], ds), file = paste(pdFiles[i], "sm", sep=""))
      }
      ## write out smoothed and non-smoothed estimates
      ds@descrip <- "R2smoothed"
      writeNIfTI(as.nifti(rcoef.s$theta[4, , , ], ds), file = "R2smoothed")
      ds@descrip <- "IPD0smoothed"
      writeNIfTI(as.nifti(rcoef.s$theta[3, , , ], ds), file = "IPD0smoothed")
      ds@descrip <- "IMT0smoothed"
      writeNIfTI(as.nifti(rcoef.s$theta[2, , , ], ds), file = "IMT0smoothed")
      ds@descrip <- "IT10smoothed"
      writeNIfTI(as.nifti(rcoef.s$theta[1, , , ], ds), file = "IT10smoothed")
      ds@descrip <- "R2"
      writeNIfTI(as.nifti(rcoef[4, , , ], ds), file = "R2")
      ds@descrip <- "IPD0"
      writeNIfTI(as.nifti(rcoef[3, , , ], ds), file = "IPD0")
      ds@descrip <- "IMT0"
      writeNIfTI(as.nifti(rcoef[2, , , ], ds), file = "IMT0")
      ds@descrip <- "IT10"
      writeNIfTI(as.nifti(rcoef[1, , , ], ds), file = "IT10")
      
    }
       
  )
    
  ## Define the value of the list within the current environment.
  assign('this', me, envir=thisEnv)
  
  ## Set the name for the class
  class(me) <- append(class(me), "mpmdata")
  return(me)
}

readMPMdata <- function(t1Files  = NULL,
                        pdFiles  = NULL,
                        mtFiles  = NULL,
                        maskFile = NULL,
                        sdim     = NULL,
                        verbose  = TRUE) {
  
  ## we need at least T1w and PDw files
  if (is.null(t1Files)) stop("vector of T1 files required")
  if (is.null(pdFiles)) stop("vector of PD files required")
  ## TODO: test whether there are enough files for the model?
  
  if (is.null(sdim)) stop("need spatial dimensionality of the data")
  if (!is.numeric(sdim) | length(sdim) != 3) stop("need exactly three numbers for spatial dimensions")
    
  ## select the model according to the existence of MTw files
  model <- if (is.null(mtFiles)) {
    1 # the simple model without MT inclusion
  } else {
    2 # the model including MT
  }
  
  ## count the number of data volumes (if (is.null(mtFiles)) length(mtFiles) == 0)
  nFiles <- length(t1Files) + length(mtFiles) + length(pdFiles)
  
  # the array for the data itself
  ddata <- numeric(nFiles * prod(sdim))
  dim(ddata) <- c(nFiles, prod(sdim))
  
  ## for each files we have a TR, TE, and flip angle (FA)
  TR <- TE <- FA <- numeric(nFiles)
  
  ## ... now we read all data volumes and extract the TR, TE, and FA values for each ...
  ii <- 1
  ## ... for all T1 volumes ...
  if (verbose) cat("reading T1 files\n")
  if (verbose) pb <- txtProgressBar(min = 0, max = length(t1Files), style = 3)
  for (i in 1:length(t1Files)) {
    ds <- readNIfTI(t1Files[i], reorient = FALSE)
    ddata[ii, ] <- ds
    ## IMPORTANT: This is special to Siawoosh data
    res <- str_match_all(ds@descrip, "([[:alpha:]]{2})=([.0123456789]+)([[:alpha:]]{2,})")
    for (nn in 1:dim(res[[1]])[1]) {
      if (res[[1]][nn, 2] == "TR") TR[ii] <- as.numeric(res[[1]][nn, 3])
      if (res[[1]][nn, 2] == "TE") TE[ii] <- as.numeric(res[[1]][nn, 3])
      if (res[[1]][nn, 2] == "FA") FA[ii] <- as.numeric(res[[1]][nn, 3])
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
      ## IMPORTANT: This is special to Siawoosh data
      res <- str_match_all(ds@descrip, "([[:alpha:]]{2})=([.0123456789]+)([[:alpha:]]{2,})")
      for (nn in 1:dim(res[[1]])[1]) {
        if (res[[1]][nn, 2] == "TR") TR[ii] <- as.numeric(res[[1]][nn, 3])
        if (res[[1]][nn, 2] == "TE") TE[ii] <- as.numeric(res[[1]][nn, 3])
        if (res[[1]][nn, 2] == "FA") FA[ii] <- as.numeric(res[[1]][nn, 3])
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
    ## IMPORTANT: This is special to Siawoosh data
    res <- str_match_all(ds@descrip, "([[:alpha:]]{2})=([.0123456789]+)([[:alpha:]]{2,})")
    for (nn in 1:dim(res[[1]])[1]) {
      if (res[[1]][nn, 2] == "TR") TR[ii] <- as.numeric(res[[1]][nn, 3])
      if (res[[1]][nn, 2] == "TE") TE[ii] <- as.numeric(res[[1]][nn, 3])
      if (res[[1]][nn, 2] == "FA") FA[ii] <- as.numeric(res[[1]][nn, 3])
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
  
  invisible(new("mpmData",
                .Data = ddata,
                t1Files = t1Files,
                mtFiles = mtFiles,
                pdFiles = pdFiles,
                TR = TR,
                TE = TE,
                FA = FA,
                smoothPar = c(0, 1, 1),
                maskFile = maskFile,
                mask = mask,
                model = as.integer(model),
                sdim = as.integer(sdim),
                nFiles = as.integer(nFiles)))
  
}

estimateESTATICS <- function(object,  ...) cat("No ESTATICS estimation defined for this class:", class(object), "\n")

setGeneric("estimateESTATICS", function(object,  ...) standardGeneric("estimateESTATICS"))

setMethod("estimateESTATICS", 
          "mpmData",
           function(object, 
                    TEScale = 100,
                    dataScale = 1000,
                    verbose = TRUE) {
             
             ## this is our design ...
             if (object@model == 2) {
               xmat <- matrix(0, object@nFiles, 4)
               xmat[1:length(object@t1Files), 1] <- 1
               xmat[(length(object@t1Files)+1):(length(object@t1Files)+length(object@mtFiles)), 2] <- 1
               xmat[(length(object@t1Files)+length(object@mtFiles)+1):(length(object@t1Files)+length(object@mtFiles)+length(object@pdFiles)), 3] <- 1
               xmat[, 4] <- object@TE / TEScale
               npar <- 4
               ## ... for our model in qflashpl() ...
               ## S_{T1} = par[1] * exp(- par[4] * TE)
               ## S_{MT} = par[2] * exp(- par[4] * TE)
               ## S_{PD} = par[3] * exp(- par[4] * TE)
             } else {
               xmat <- matrix(0, object@nFiles, 3)
               xmat[1:length(object@t1Files), 1] <- 1
               xmat[(length(object@t1Files)+1):(length(object@t1Files)+length(object@pdFiles)), 2] <- 1
               xmat[, 3] <- object@TE / TEScale
               npar <- 3
               ## ... for our model in qflashpl2() ...
               ## S_{T1} = par[1] * exp(- par[3] * TE)
               ## S_{PD} = par[2] * exp(- par[3] * TE)        
             }
             if (verbose) {
               cat("Design of the model:\n")
               print(xmat)
             }
             
             ## starting value for R* estimate
             R2star <- 0.05 * TEScale
             indT1 <- order(object@TE[as.logical(xmat[, 1])])[1]
             if (object@model == 2) {
               indMT <- order(object@TE[as.logical(xmat[, 2])])[1] + sum(xmat[, 1])
               indPD <- order(object@TE[as.logical(xmat[, 3])])[1] + sum(xmat[, 1]) + sum(xmat[, 2]) 
             } else {
               indPD <- order(object@TE[as.logical(xmat[, 2])])[1] + sum(xmat[, 1])        
             }
  
             modelCoeff <- new("ESTATICSModel",
                               .Data = array(0, c(npar, object@sdim)),
                               isConv = array(FALSE, object@sdim),
                               invCov = array(0, c(npar, npar, object@sdim)),
                               bi = array(1, object@sdim),
                               TEScale = TEScale,
                               dataScale = dataScale,
                               smoothPar = c(0, 1, 1),
                               mask = object@mask,
                               model = object@model,
                               sdim = object@sdim,
                               nFiles = object@nFiles)
             
             ## now perform the voxelwise regression
             if (verbose) Sys.time()
             for (z in 1:object@sdim[3]) {
               for (y in 1:object@sdim[2]) {
                 for (x in 1:object@sdim[1]) {
                   if (object@mask[x, y, z]) {
                     ivec  <- object[, x, y, z] / dataScale
                     if (object@model == 2) { 
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
                     modelCoeff@isConv[x, y, z] <- res$convInfo$isConv
                     modelCoeff@.Data[, x, y, z] <- sres$coefficients[, 1]
                     if (sres$sigma != 0) {
                       modelCoeff@invCov[, , x, y, z] <- solve(sres$cov.unscaled) / sres$sigma^2
                     } else {
                       ## TODO: CHECK!
                       modelCoeff@invCov[, , x, y, z] <- 0
                     }
                   }
                 }
               }
               if (verbose) cat(z, format(Sys.time()), "\n")
             }
             if (verbose) Sys.time()
             
             ## END function estimateESTATICS()
             return(invisible(modelCoeff))
           })

smoothESTATICS <- function(object,  ...) cat("No ESTATICS smoothing defined for this class:", class(object), "\n")

setGeneric("smoothESTATICS", function(object,  ...) standardGeneric("smoothESTATICS"))

setMethod("smoothESTATICS", 
          "ESTATICSModel",
          function(object, 
                   mpmData = NULL,
                   smoothedDataFile = NULL,
                   kstar = 16, 
                   alpha = 0.05, 
                   wghts = NULL,
                   verbose = TRUE) {

            ## length of the vector to smooth (# parameters of model)
            nv <- if (object@model == 2) 4 else 3
            
            ## determine a suitable adaptation bandwidth
            lambda <- nv * qf(1 - alpha, nv, object@nFiles - nv)
            
            ## adjust for non-isotropic voxel if necessary
            if (is.null(wghts)) wghts <- c(1, 1, 1)
            ## make first spatial dimension unit, and use second and third for reference
            wghts <- wghts[1] / wghts[2:3]
            
            ## spatial dimension and number of voxel
            n1 <- object@sdim[1]
            n2 <- object@sdim[2]
            n3 <- object@sdim[3]
            n <- n1 * n2 * n3
            
            ## initialization for first step
            zobj <- list(bi = rep(1, n), theta = object@.Data)
            bi <- zobj$bi
            
            ## if we plan to smooth the original data too, we take special care in the last step
            if (!is.null(mpmData)) {
              if (is.null(smoothedDataFile)) stop("need file name for smoothed data")
              if (file.exists(smoothedDataFile)) {
                cont <- readline("file", smoothedDataFile, "already exists. Continue [Y/N]?")
                if (toupper(cont) != "Y") stop("stopping")
              }
              if(length(mpmData) != 4 | any(object@sdim != dim(mpmData)[-1]))  
                stop("incompatible dimensions of model parameters and original data")
              kstar <- kstar-1
            }
            
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
            while (k <= kstar) {
              ## determine the actual bandwidth for this step
              hakt <- gethani(1, 1.25*hmax, 2, 1.25^k, wghts, 1e-4)
              
              ## we need the (approx.) size of the weigthing scheme array 
              dlw <- (2*trunc(hakt/c(1, wghts))+1)[1:3]
              
              ## perform the actual adaptive smoothing 
              zobj <- .Fortran("vaws",
                               as.double(object@.Data),
                               as.logical(object@mask),
                               as.integer(nv),
                               as.integer(n1),
                               as.integer(n2),
                               as.integer(n3),
                               hakt = as.double(hakt),
                               as.double(lambda),
                               as.double(zobj$theta),
                               as.double(object@invCov),
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
                                     "MSE: ", signif(mean((zobj$theta - object@.Data)^2),3),
                                     "MAE: ", m1 <- signif(mean(abs(zobj$theta - object@.Data)),3),
                                     "mean(bi):", signif(mean(zobj$bi),3),
                                     "\n")
                mae <- c(mae, m1)
                setTxtProgressBar(pb, i)
              }
              
              ## go for next iteration
              k <- k+1
              gc()
              
            }
            if (!is.null(mpmData)) { ##  modified last step; smoothing data, too
              
              ## we need the number of files for some array dimensions
              nve <- object@nFiles
              
              ## determine the actual bandwidth for this step
              hakt <- gethani(1, 1.25*hmax, 2, 1.25^k, wghts, 1e-4)
              
              ## we need the (approx.) size of the weigthing scheme array 
              dlw <- (2*trunc(hakt/c(1, wghts))+1)[1:3]
              
              ## perform the actual adaptive smoothing
              zobj <- .Fortran("vawsext",
                               as.double(object@.Data),
                               as.logical(object@mask),
                               as.integer(nv),
                               as.integer(n1),
                               as.integer(n2),
                               as.integer(n3),
                               as.double(mpmData@.Data),
                               as.integer(nve),
                               hakt = as.double(hakt),
                               as.double(lambda),
                               as.double(zobj$theta),
                               as.double(object@invCov),
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
              dim(zobj$thext) <- c(object@nFiles, object@sdim)
              
              mpmDatasmoothed <- new("mpmData",
                                     .Data = zobj$thext,
                                     t1Files = mpmData@t1Files,
                                     mtFiles = mpmData@mtFiles,
                                     pdFiles = mpmData@pdFiles,
                                     TR = mpmData@TR,
                                     TE = mpmData@TE,
                                     FA = mpmData@FA,
                                     smoothPar = c(lambda, hakt, alpha),
                                     maskFile = mpmData@maskFile,
                                     mask = mpmData@mask,
                                     model = mpmData@model,
                                     sdim = mpmData@sdim,
                                     nFiles = mpmData@nFiles)
              save(mpmDatasmoothed, file = smootheDataFile)
              
              ## use maximum ni
              zobj$bi <- pmax(bi, zobj$bi)
              
              ## some verbose stuff
              if (verbose) {
                protocol[k] <- paste("bandwidth: ", signif(hakt, 3),
                                     "MSE: ", signif(mean((zobj$theta - get("modelCoeff", thisEnv))^2),3),
                                     "MAE: ", m1 <- signif(mean(abs(zobj$theta - get("modelCoeff", thisEnv))),3),
                                     "mean(bi):", signif(mean(zobj$bi),3),
                                     "\n")
                mae <- c(mae, m1)
                setTxtProgressBar(pb, i)
              }
                            
            }
            if (verbose) close(pb)
            if (verbose) print(protocol)
            
            ## set correct dimensions
            dim(zobj$theta) <- c(nv, n1, n2, n3)
            dim(zobj$bi) <- c(n1, n2, n3)
            
            ## assign values
            modelCoeffsmooth <- new("ESTATICSModel",
                                    .Data = zobj$theta,
                                    isConv = object@isConv,
                                    invCov = object@invCov,
                                    bi = zobj$bi,
                                    TEScale = object@TEScale,
                                    dataScale = object@dataScale,
                                    smoothPar = c(lambda, hakt, alpha),
                                    mask = object@mask,
                                    model = object@model,
                                    sdim = object@sdim,
                                    nFiles = object@nFiles)
            
            return(invisible(modelCoeffsmooth))
            ## END function smoothESTATICS()  
          })



                      