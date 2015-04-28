mpmData <- function(myt1Files  = NULL,
                    mypdFiles  = NULL,
                    mymtFiles  = NULL,
                    mymaskFile = NULL,
                    myb1File   = NULL,
                    mysdim       = NULL) {
  
  
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
    1 # the simple model with MT inclusion ("second cosine")
  } else {
    2 # the model including MT
  }
  
  ## count the number of data volumes (if (is.null(mtFiles)) length(mtFiles) == 0)
  nFiles <- length(t1Files) + length(mtFiles) + length(pdFiles)

  sdim <- mysdim
  
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
    
    getB1File = function()
    {
      return(get("b1File", thisEnv))
    },
    
    readMPMData = function(verbose = FALSE) {

      ## for each files we have a TR, TE, and flip angle (FA)
      t1TR <- t1TE <- t1FA <- numeric(length(get("t1Files", thisEnv)))
      pdTR <- pdTE <- pdFA <- numeric(length(get("pdFiles", thisEnv)))
      if (get("model", thisEnv) == 2) {
        mtTR <- mtTE <- mtFA <- numeric(length(get("mtFiles", thisEnv)))
      }

      ## we want to read all data into a 4 dimensional data cube ...
      if (get("model", thisEnv) == 2) {
        ddim <- c(length(t1TR) + length(mtTR) + length(pdTR), get("sdim", thisEnv))
      } else {
        ddim <- c(length(t1TR) + length(pdTR), get("sdim", thisEnv))        
      }
      dd <- numeric(prod(ddim))
      dim(dd) <- c(ddim[1], prod(ddim[2:4]))
      
      ## ... now we read all data volumes and extract the TR, TE, and FA values for each ...
      ii <- 1
      ## ... for all T1 volumes ...
      if (verbose) cat("reading T1 files\n")
      if (verbose) pb <- txtProgressBar(min = 0, max = length(get("t1Files", thisEnv)), style = 3)
      for (i in 1:length(get("t1Files", thisEnv))) {
        ds <- readNIfTI(get("t1Files", thisEnv)[i], reorient=FALSE)
        dd[ii, ] <- ds
        res <- str_match_all(ds@descrip, "([[:alpha:]]{2})=([.0123456789]+)([[:alpha:]]{2,})")
        for (nn in 1:dim(res[[1]])[1]) {
          if (res[[1]][nn, 2] == "TR") t1TR[i] <- as.numeric(res[[1]][nn, 3])
          if (res[[1]][nn, 2] == "TE") t1TE[i] <- as.numeric(res[[1]][nn, 3])
          if (res[[1]][nn, 2] == "FA") t1FA[i] <- as.numeric(res[[1]][nn, 3])
        }
        ii <- ii + 1
        if (verbose) setTxtProgressBar(pb, i)
      }
      if (verbose) close(pb)
      ## ... for all MT volumes ...
      if (get("model", thisEnv) == 2) {
        if (verbose) cat("reading MT files\n")
        if (verbose) pb <- txtProgressBar(min = 0, max = length(get("mtFiles", thisEnv)), style = 3)
        for (i in 1:length(get("mtFiles", thisEnv))) {
          ds <- readNIfTI(get("mtFiles", thisEnv)[i], reorient=FALSE)
          dd[ii, ] <- ds
          res <- str_match_all(ds@descrip, "([[:alpha:]]{2})=([.0123456789]+)([[:alpha:]]{2,})")
          for (nn in 1:dim(res[[1]])[1]) {
            if (res[[1]][nn, 2] == "TR") mtTR[i] <- as.numeric(res[[1]][nn, 3])
            if (res[[1]][nn, 2] == "TE") mtTE[i] <- as.numeric(res[[1]][nn, 3])
            if (res[[1]][nn, 2] == "FA") mtFA[i] <- as.numeric(res[[1]][nn, 3])
          }
          ii <- ii + 1
          if (verbose) setTxtProgressBar(pb, i)
        }
      }
      if (verbose) close(pb)
      ## .. and for all PD volumes ...
      if (verbose) cat("reading PD files\n")
      if (verbose) pb <- txtProgressBar(min = 0, max = length(get("pdFiles", thisEnv)), style = 3)
      for (i in 1:length(get("pdFiles", thisEnv))) {
        ds <- readNIfTI(get("pdFiles", thisEnv)[i], reorient=FALSE)
        dd[ii, ] <- ds
        res <- str_match_all(ds@descrip, "([[:alpha:]]{2})=([.0123456789]+)([[:alpha:]]{2,})")
        for (nn in 1:dim(res[[1]])[1]) {
          if (res[[1]][nn, 2] == "TR") pdTR[i] <- as.numeric(res[[1]][nn, 3])
          if (res[[1]][nn, 2] == "TE") pdTE[i] <- as.numeric(res[[1]][nn, 3])
          if (res[[1]][nn, 2] == "FA") pdFA[i] <- as.numeric(res[[1]][nn, 3])
        }
        ii <- ii + 1
        if (verbose) setTxtProgressBar(pb, i)
      }
      if (verbose) close(pb)
      dim(dd) <- ddim
      ## ... done!
      
      ## TODO: Does this involve some unneccessary copy of a large array?
      ## TODO: -> create it first and then assign only parts.
      assign("dd", dd, thisEnv)
      if (get("model", thisEnv) == 2) {
        assign("TR", c(t1TR, mtTR, pdTR), thisEnv)
      } else {
        assign("TR", c(t1TR, pdTR), thisEnv)
      }
      if (get("model", thisEnv) == 2) {
        assign("TE", c(t1TE, mtTE, pdTE), thisEnv)
      } else {
        assign("TE", c(t1TE, pdTE), thisEnv)
      }
      if (get("model", thisEnv) == 2) {
        assign("FA", c(t1FA, mtFA, pdFA), thisEnv)
      } else {
        assign("FA", c(t1FA, pdFA), thisEnv)
      }
    },
    
    getMPMData = function() {
      return(get("dd", thisEnv))
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
    
    estimateESTATICS = function() {
      ## we need some scales to imrpove the condition number of the covariance matrix:
      TEscale <- 100
      ivecScale <- 1000
      
      ## this is our design ...
      xmat <- matrix(0, ddim[1], 4)
      xmat[1:length(t1Files), 1] <- 1
      xmat[(length(t1Files)+1):(length(t1Files)+length(mtFiles)), 2] <- 1
      xmat[(length(t1Files)+length(mtFiles)+1):(length(t1Files)+length(mtFiles)+length(pdFiles)), 3] <- 1
      xmat[, 4] <- c(t1TE, mtTE, pdTE) / TEscale
      ## ... for our model in qflashpl() ...
      ## I_{T1} = par[1] * exp(- par[4] * TE)
      ## I_{MT} = par[2] * exp(- par[4] * TE)
      ## I_{PD} = par[3] * exp(- par[4] * TE)
      library(qMRI)
      
      
      ## read brain mask file (created by FSL bet2 using 0.3 as fractional intensity threshold)
      mask <- as.logical(readNIfTI(maskFile, reorient = FALSE))
      dim(mask) <- c(280, 320, 208)
      
      ## now perform the voxelwise regression
      Sys.time()
      rcoef <- array(0, c(4, 280, 320, 208))
      rrcov <- array(0, c(4, 4, 280, 320, 208))
      rsigma <- array(0, c(280, 320, 208))
      isConv <- array(FALSE, c(280, 320, 208))
      pNames <- c("ST1","SMT","SPD","R2star")
      dimnames(rcoef) <- list(pNames,NULL,NULL,NULL)
      dimnames(rrcov) <- list(pNames,pNames,NULL,NULL,NULL)
      for (z in 1:208){
        for (y in 1:320) {
          for (x in 1:280) {
            if (mask[x, y, z]) {
              ivec  <- dd[, x, y, z] / ivecScale
              th <- c(ivec[1]*exp(-xmat[1, 4]*5),   # par[1]
                      ivec[9]*exp(-xmat[9, 4]*5),   # par[2]
                      ivec[15]*exp(-xmat[15, 4]*5), # par[3]
                      5)                            # par[4]
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
              sres <- summary(res) 
              isConv[x, y, z] <- res$convInfo$isConv
              rrcov[, , x, y, z] <- sres$cov.unscaled
              rsigma[x, y, z] <- sres$sigma
              rcoef[, x, y, z] <- sres$coefficients[, 1]
            }
          }
        }
        cat(z, format(Sys.time()), "\n")
        #  save(rrcov, rcoef, rsigma, isConv, file="mpmESTATICS.rsc")
      }
      Sys.time()
      save(rrcov, rcoef, rsigma, isConv, file="mpmESTATICS.rsc")
      
      
      ## we now determine the inverse covariance matrix for adaptive smoothing
      load("mpmESTATICS.rsc")
      dim(rrcov) <- c(4, 4, 280*320*208)
      irrcov <- array(0, dim(rrcov))
      irrcov[, , mask] <- apply(rrcov[, , mask], 3, solve)
      dim(rsigma) <- NULL
      irrcovF <- sweep(irrcov, 3, rsigma^2, "/")
      dim(irrcovF) <- dim(irrcov) <- c(4, 4, 280, 320, 208)
      irrcovF[is.na(irrcovF)] <- 0
      dim(rrcov) <- c(4, 4, 280, 320, 208)
      save(irrcov, irrcovF, file="irrcovESTATICS.rsc")
      
    },
    
    smoothESTATICS = function() {
      library(qMRI)
      setCores(10)
      rcoef.s <- vaws(rcoef, 
                      irrcovF, 
                      kstar = 16, 
                      mask = mask, 
                      yext = dd,
                      alpha = 0.05, 
                      df = 18, 
                      ladjust = 1, 
                      wghts = NULL, 
                      u = rcoef, 
                      maxni = TRUE)
      
    },
    
    getQMapsFromESTATICS <- function(mpmData, b1Map) {
      
      b1Map <- readNIfTI(b1File, reorient = FALSE)/100
      b1Map[b1Map < 0] <- 0
      
      ## calculate f1
      alphat1 <- b1Map * t1FA[1] /180 * pi
      alphapd <- b1Map * pdFA[1] /180 * pi
      SINalphat1 <- sin(alphat1)
      COSalphat1 <- cos(alphat1)
      SINalphapd <- sin(alphapd)
      COSalphapd <- cos(alphapd)
      rm(alphat1, alphapd)
      enum <- rcoef.s$theta[1, , , ] - SINalphat1/SINalphapd * rcoef.s$theta[3, , , ]
      denom <- rcoef.s$theta[1, , , ] * COSalphat1 - SINalphat1/SINalphapd * rcoef.s$theta[3, , , ] * COSalphapd
      f1 <- enum/denom
      rm(enum, denom, COSalphapd, SINalphapd)
      R1 <- -log(f1)/t1TR[1]
      
      ## calculate PD
      enum <- (1 - COSalphat1 * f1) * rcoef.s$theta[1, , , ]
      denom <- SINalphat1 * (1 - f1)
      PD <- enum/denom
      rm(enum, denom, SINalphat1)
      
      ## calculate delta
      alphamt <- b1Map * mtFA[1] /180 * pi
      denom1 <- rcoef.s$theta[2, , , ] * f1 * cos(alphamt)
      denom2 <- PD * (1 - f1) * sin(alphamt)
      delta <- 1 - rcoef.s$theta[4, , , ]/(denom1 + denom2)
      rm(alphamt, denom1, denom2)
      
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



