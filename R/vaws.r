gethani <- function(x,y,lkern,value,wght,eps=1e-2){
  .Fortran(C_gethani,
           as.double(x),
           as.double(y),
           as.integer(lkern),
           as.double(value),
           as.double(wght),
           as.double(eps),
           bw=double(1))$bw
}

geth.gauss <- function(corr, step = 1.01) {
  #   get the   bandwidth for lkern corresponding to a given correlation
  #
  #  keep it simple result does not depend on d
  #
  if (corr < 0.1) {
    h <- 1e-5
  } else {
    h <- .8
    z <- 0
    while (z < corr) {
      h <- h * step
      z <- get.corr.gauss(h, interv = 2)
    }
    h <- h / step
  }
  h
}
get.corr.gauss <- function(h, interv = 1) {
  #
  #   Calculates the correlation of
  #   colored noise that was produced by smoothing with "gaussian" kernel and bandwidth h
  #   Result does not depend on d for "Gaussian" kernel !!
  h <- h / 2.3548 * interv
  ih <- trunc(4 * h + 2 * interv - 1)
  dx <- 2 * ih + 1
  penl <- dnorm(((-ih):ih) / h)
  sum(penl[-(1:interv)] * penl[-((dx - interv + 1):dx)]) / sum(penl ^
                                                                 2)
}

smoothInvCov <- function(mpmESTATICSobj, hcov=3, verbose=TRUE,
   method=c("rss","logInvCov","InvCov")){
    sdim <- mpmESTATICSobj$sdim
    mask <- mpmESTATICSobj$mask
    if(method=="rss"){
       rss <- mpmESTATICSobj$rsigma^2
       rss[!mask] <- mean(rss[mask])
## avoid effects from 0 values outside mask
       rsssm <- kernsm(rss, hcov)@yhat
       rsssm[rsssm<=0] <- rss[rsssm<=0]
       qrss <- rss/rsssm
       rsssm[!mask] <- 0
       qrss[!mask] <- 1
       mpmESTATICSobj$invCov <- sweep(mpmESTATICSobj$invCov,3:5,qrss,"*")
       mpmESTATICSobj$rsigma <- sqrt(rsssm)
    }
    if(method=="logInvCov"){
      logInvCov <- mpmESTATICSobj$invCov
      d <- dim(logInvCov)[1]
      for (z in 1:sdim[3]){
        for (y in 1:sdim[2]) {
          for (x in 1:sdim[1]) {
            if (mask[x, y, z]) {
               zsvd <- svd(logInvCov[,,x,y,z])
               zd <- zsvd$d
               zd <- pmax(zd,1e-8*zd[1])
               logInvCov[,,x,y,z] <- zsvd$u%*%diag(log(zd))%*%t(zsvd$u)
            }
          }
        }
      if(verbose) cat(".")
      }
      if(verbose) cat("\n start smoothing \n")
      for(i in 1:d){
        for(j in 1:i){
           covij <- logInvCov[i,j,,,]
           covij[!mask] <- mean(covij[mask])
           ## avoid effects from 0 values outside mask
           zsm <- kernsm(covij, hcov)@yhat
           zsm[!mask] <- 0
           logInvCov[i,j,,,] <- zsm
           if(i>j) logInvCov[j,i,,,] <- zsm
        }
      }
      for (z in 1:sdim[3]){
        for (y in 1:sdim[2]) {
          for (x in 1:sdim[1]) {
            if (mask[x, y, z]) {
               zsvd <- svd(logInvCov[,,x,y,z])
               logInvCov[,,x,y,z] <- zsvd$u%*%diag(exp(zsvd$d))%*%t(zsvd$u)
            }
          }
        }
        if(verbose) cat(".")
      }
      mpmESTATICSobj$invCov <- logInvCov
    }
  if(method=="InvCov"){
     InvCov <- mpmESTATICSobj$invCov
     d <- dim(InvCov)[1]
     if(verbose) cat("\n start smoothing \n")
     for(i in 1:d){
       for(j in 1:i){
         covij <- InvCov[i,j,,,]
         covij[!mask] <- mean(covij[mask])
         ## avoid effects from 0 values outside mask
          zsm <- kernsm(covij, hcov)@yhat
          zsm[!mask] <- 0
          InvCov[i,j,,,] <- zsm
          if(i>j) InvCov[j,i,,,] <- zsm
       }
     }
     mpmESTATICSobj$invCov <- InvCov
  }
    mpmESTATICSobj
}


vpawscov0 <- function(y,
                     kstar = 16,
                     invcov = NULL,
                     mask = NULL,
                     scorr = 0,
                     spmin = 0.25,
                     lambda = NULL,
                     ladjust = 1,
                     wghts = NULL,
                     maxni = FALSE,
                     patchsize = 1,
                     data = NULL,
                     verbose = TRUE) {
  ##
  ##  this is the version with full size invcov (triangular storage)
  ##  needed for MPM
  ##  version for complete 3D cubes
  ##
  dy <- dim(y)
  nvec <- dy[1]
  if(nvec>5) stop("limited to 5 parameters")
  indcov <- switch(nvec,1,
                   c(1,2,4),
                   c(1,2,5,3,6,9),
                   c(1,2,6,3,7,11,4,8,12,16),
                   c(1,2,7,3,8,13,4,9,14,19,5,10,15,20,25))
  if(!is.null(data)) nsample <- dim(data)[1]
  dy <- dy[-1]
  d <- length(dy)
  if (d != 3)
    stop("need 3 dimensional grids")
  if(is.null(lambda)){
    lambda <- 2 * ladjust * qchisq(pchisq(8.82, 1), nvec)
  }
  if (is.null(wghts)) wghts <- c(1, 1, 1)
  wghts <- wghts[1] / wghts[2:3]
  n1 <- dy[1]
  n2 <- dy[2]
  n3 <- dy[3]
  n <- n1 * n2 * n3
  if (is.null(mask))
    mask <- rep(TRUE, n)
  h0 <- 0
  if (any(scorr > 0)) {
    h0 <- numeric(length(scorr))
    for (i in 1:length(h0))
      h0[i] <- geth.gauss(scorr[i])
    if (length(h0) < d)
      h0 <- rep(h0[1], d)
    if(verbose) cat("Corresponding bandwiths for specified correlation:",
        h0,
        "\n")
  }
  ## create index information for voxel in mask
  nvoxel <- sum(mask)
  position <- array(0,dy)
  position[mask] <- 1:nvoxel
  dim(mask) <- NULL
  dim(y) <- c(nvec,n)
  dim(invcov) <- c(nvec * nvec,n)
  hseq <- 1
  zobj <- list(bi = rep(1, nvoxel), theta = y[,as.vector(mask)])
  bi <- zobj$bi
  if(verbose) cat("Progress:")
  total <- cumsum(1.25 ^ (1:kstar)) / sum(1.25 ^ (1:kstar))
  mc.cores <- setCores(, reprt = FALSE)
  np1 <- 2 * patchsize + 1
  np2 <- if (n2 > 1) 2 * patchsize + 1 else 1
  np3 <- if (n3 > 1) 2 * patchsize + 1 else 1
  k <- 1
  hmax <- 1.25 ^ (kstar / d)
  lambda0 <- 1e32
  mae <- NULL
  while (k <= kstar) {
    hakt0 <- gethani(1, 1.25 * hmax, 2, 1.25 ^ (k - 1), wghts, 1e-4)
    hakt <- gethani(1, 1.25 * hmax, 2, 1.25 ^ k, wghts, 1e-4)
    if(verbose) cat("step", k, "hakt", hakt, "time", format(Sys.time()), "\n")
    hseq <- c(hseq, hakt)
    dlw <- (2 * trunc(hakt / c(1, wghts)) + 1)
    if(k==kstar & !is.null(data)){
      dim(data) <- c(nsample,n)
      zobj <- .Fortran(C_pvawsme,
                       as.double(y[,mask]),
                       as.double(data[,mask]), ## data to smooth additionally
                       as.integer(position),
                       as.integer(nvec),
                       as.integer(nvec * (nvec + 1) / 2),
                       as.integer(nsample), ## leading dimension of data
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt = as.double(hakt),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       as.double(zobj$bi),
                       bi = double(nvoxel), #binn
                       theta = double(nvec * nvoxel),
                       data = double(nsample*nvoxel),
                       as.double(invcov[indcov,mask]),#
                       as.integer(mc.cores),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nvec * mc.cores),
                       double(nsample * mc.cores),
                       as.integer(np1),
                       as.integer(np2),
                       as.integer(np3))[c("bi", "theta", "hakt","data")]
      data[,mask] <- zobj$data
      dim(data) <- c(nsample, dy)
    } else {
      zobj <- .Fortran(C_pvawsm2,
                       as.double(y[,mask]),
                       as.integer(position),
                       as.integer(nvec),
                       as.integer(nvec * (nvec + 1) / 2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt = as.double(hakt),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       as.double(zobj$bi),
                       bi = double(nvoxel), #binn
                       theta = double(nvec * nvoxel),
                       as.double(invcov[indcov,mask]),# compact storage
                       as.integer(mc.cores),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nvec * mc.cores),
                       as.integer(np1),
                       as.integer(np2),
                       as.integer(np3))[c("bi", "theta", "hakt")]
    }
    if (maxni)
      bi <- zobj$bi <- pmax(bi, zobj$bi)
    x <- 1.25 ^ k
    scorrfactor <- x / (3 ^ d * prod(scorr) * prod(h0) + x)
    lambda0 <- lambda * scorrfactor
    if (verbose & max(total) > 0) {
      cat(signif(total[k], 2) * 100, "%  ", sep = "")
      cat("mean(bi)", signif(mean(zobj$bi),3)," ")
    }
    k <- k + 1
    gc()
  }
  theta <- matrix(0,nvec,n)
  theta[,mask] <- zobj$theta
  dim(theta) <- c(nvec, dy)
  bi <- numeric(n)
  bi[mask] <- zobj$bi
  dim(bi) <- dy
  if(verbose) cat("\n")
  list(
    theta=theta,
    hakt=hakt,
    lambda=lambda,
    hseq = hseq,
    bi = bi,
    data= if(!is.null(data)) data else NULL
  )
}

vpawscov <- function(y,
                     kstar = 16,
                     invcov = NULL,
                     mask = NULL,
                     scorr = 0,
                     spmin = 0.25,
                     lambda = NULL,
                     ladjust = 1,
                     wghts = NULL,
                     maxni = FALSE,
                     patchsize = 1,
                     data = NULL,
                     verbose = TRUE) {#1
  ##
  ##  this is the version with full size invcov (triangular storage)
  ##  needed for MPM
  ##  version for data reduced to mask
  ##
  dy <- dim(y)
  nvec <- dy[1]
  if(nvec>5) stop("limited to 5 parameters")
  indcov <- switch(nvec,1,
                   c(1,2,4),
                   c(1,2,5,3,6,9),
                   c(1,2,6,3,7,11,4,8,12,16),
                   c(1,2,7,3,8,13,4,9,14,19,5,10,15,20,25))
  if(!is.null(data)) nsample <- dim(data)[1]
  dy <- dim(mask)
  d <- length(dy)
  if (d != 3)
    stop("need 3D mask")
  if(is.null(lambda)){#2
    lambda <- 2 * ladjust * qchisq(pchisq(8.82, 1), nvec)
  }#2
  if (is.null(wghts)) wghts <- c(1, 1, 1)
  wghts <- wghts[1] / wghts[2:3]
  n1 <- dy[1]
  n2 <- dy[2]
  n3 <- dy[3]
  h0 <- 0
  if (any(scorr > 0)) {#3
    h0 <- numeric(length(scorr))
    for (i in 1:length(h0))
      h0[i] <- geth.gauss(scorr[i])
    if (length(h0) < d)
      h0 <- rep(h0[1], d)
    if(verbose) cat("Corresponding bandwiths for specified correlation:",
        h0,
        "\n")
  }#3
  ## create index information for voxel in mask
  nvoxel <- sum(mask)
  position <- array(0,dy)
  position[mask] <- 1:nvoxel
  dim(mask) <- NULL
  dim(y) <- c(nvec,nvoxel)
  dim(invcov) <- c(nvec * nvec,nvoxel)
  hseq <- 1
  zobj <- list(bi = rep(1, nvoxel), theta = y)
  bi <- zobj$bi
  if(verbose) cat("Progress:")
  total <- cumsum(1.25 ^ (1:kstar)) / sum(1.25 ^ (1:kstar))
  mc.cores <- setCores(, reprt = FALSE)
  np1 <- 2 * patchsize + 1
  np2 <- if (n2 > 1) 2 * patchsize + 1 else 1
  np3 <- if (n3 > 1) 2 * patchsize + 1 else 1
  k <- 1
  hmax <- 1.25 ^ (kstar / d)
  lambda0 <- 1e32
  mae <- NULL
  while (k <= kstar) {#4
    hakt0 <- gethani(1, 1.25 * hmax, 2, 1.25 ^ (k - 1), wghts, 1e-4)
    hakt <- gethani(1, 1.25 * hmax, 2, 1.25 ^ k, wghts, 1e-4)
    if(verbose) cat("step", k, "hakt", hakt, "time", format(Sys.time()), "\n")
    hseq <- c(hseq, hakt)
    dlw <- (2 * trunc(hakt / c(1, wghts)) + 1)
    if(k==kstar & !is.null(data)){#5
      dim(data) <- c(nsample,nvoxel)
      zobj <- .Fortran(C_pvawsme,
                       as.double(y),
                       as.double(data), ## data to smooth additionally
                       as.integer(position),
                       as.integer(nvec),
                       as.integer(nvec * (nvec + 1) / 2),
                       as.integer(nsample), ## leading dimension of data
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt = as.double(hakt),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       as.double(zobj$bi),
                       bi = double(nvoxel), #binn
                       theta = double(nvec * nvoxel),
                       data = double(nsample*nvoxel),
                       as.double(invcov[indcov,]),#
                       as.integer(mc.cores),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nvec * mc.cores),
                       double(nsample * mc.cores),
                       as.integer(np1),
                       as.integer(np2),
                       as.integer(np3))[c("bi", "theta", "hakt","data")]
      dim(zobj$data) <- c(nsample, nvoxel)
    } else {#6
      zobj <- .Fortran(C_pvawsm2,
                       as.double(y),
                       as.integer(position),
                       as.integer(nvec),
                       as.integer(nvec * (nvec + 1) / 2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt = as.double(hakt),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       as.double(zobj$bi),
                       bi = double(nvoxel), #binn
                       theta = double(nvec * nvoxel),
                       as.double(invcov[indcov,]),# compact storage
                       as.integer(mc.cores),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nvec * mc.cores),
                       as.integer(np1),
                       as.integer(np2),
                       as.integer(np3))[c("bi", "theta", "hakt")]
    }#6
    if (maxni)
      bi <- zobj$bi <- pmax(bi, zobj$bi)
    x <- 1.25 ^ k
    scorrfactor <- x / (3 ^ d * prod(scorr) * prod(h0) + x)
    lambda0 <- lambda * scorrfactor
    if (verbose & max(total) > 0) {#7
      cat(signif(total[k], 2) * 100, "%  ", sep = "")
      cat("mean(bi)", signif(mean(zobj$bi),3)," ")
    }#7
    k <- k + 1
    gc()
  }
  dim(zobj$theta) <- c(nvec, nvoxel)
  if(verbose) cat("\n")
  list(
    theta=zobj$theta,
    hakt=hakt,
    lambda=lambda,
    hseq = hseq,
    bi = zobj$bi,
    data= if(!is.null(zobj$data)) zobj$data else NULL
  )
}
