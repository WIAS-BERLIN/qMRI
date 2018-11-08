vaws <- function(y,si2,kstar=16,mask=NULL,yext=NULL,alpha=0.001,df=18,
                 ladjust=1,wghts=NULL,u=NULL,maxni=FALSE){
#
#  vectorized aws for 3D images
#
#  y - image data, dimension c(nv,n1,n2,n3)
#  si2 - local inverse covariance matrices  dim(nv,nv,n1,n2,n3)
#
  args <- match.call()
  dy <- dim(y)
  dsi <- dim(si2)
  if(length(dy)!=4|length(dsi)!=5|any(dy!=dsi[-1])|dsi[1]!=dy[1])
      stop("incompatible dimensions of y and si2")
  nv <- dy[1]
  dy <- dy[-1]
  if(!is.null(yext)) {
    dyext <- dim(yext)
    if(length(dyext)!=4|any(dy!=dyext[-1]))
      stop("incompatible dimensions of y and yext")
  }
  lambda <- nv*qf(1-alpha,nv,df)
  if(is.null(wghts)) wghts <- c(1,1,1)
  wghts <- wghts[1]/wghts[2:3]
  n1 <- dy[1]
  n2 <- dy[2]
  n3 <- dy[3]
  n <- n1*n2*n3
  if(is.null(mask)) mask <- rep(TRUE,n)
  h0 <- 0
  zobj<-list(bi= rep(1,n), theta= y)
  bi <- zobj$bi
  cat("Progress:")
  total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
  mc.cores <- setCores(,reprt=FALSE)
  k <- 1
  hmax <- 1.25^(kstar/3)
  lambda0 <- lambda
  mae <- NULL
  if(!is.null(yext)) kstar <- kstar-1
  while (k<=kstar) {
    hakt0 <- gethani(1,1.25*hmax,2,1.25^(k-1),wghts,1e-4)
    hakt <- gethani(1,1.25*hmax,2,1.25^k,wghts,1e-4)
    cat("step",k,"hakt",hakt,"time",format(Sys.time()),"\n")
    dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:3]
    zobj <- .Fortran(C_vaws,as.double(y),
                     as.logical(mask),
                     as.integer(nv),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     hakt=as.double(hakt),
                     as.double(lambda0),
                     as.double(zobj$theta),
                     as.double(si2),
                     bi=as.double(zobj$bi),
                     theta=double(nv*n),
                     as.integer(mc.cores),
                     double(prod(dlw)),
                     as.double(wghts),
                     double(nv*mc.cores))[c("bi","theta","hakt")]
    dim(zobj$theta)<-c(nv,dy)
    if(maxni) bi <- zobj$bi <- pmax(bi,zobj$bi)
    dim(zobj$bi)<-dy
    if(!is.null(u)) {
      cat("bandwidth: ",signif(hakt,3),"   MSE: ",
          signif(mean((zobj$theta-u)^2),3),"   MAE: ",
          signif(mean(abs(zobj$theta-u)),3)," mean(bi)=",
          signif(mean(zobj$bi),3),"\n")
      mae<-c(mae,signif(mean(abs(zobj$theta-u)),3))
    }
    if (max(total) >0) {
      cat(signif(total[k],2)*100,"% . ",sep="")
    }
    k <- k+1
    gc()
  }
  if(!is.null(yext)){
##
##  modified last step smoothing also yext
##
    nve <- dyext[1]
    hakt0 <- gethani(1,1.25*hmax,2,1.25^(k-1),wghts,1e-4)
    hakt <- gethani(1,1.25*hmax,2,1.25^k,wghts,1e-4)
    cat("step",k,"hakt",hakt,"time",format(Sys.time()),"\n")
    dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:3]
    zobj <- .Fortran(C_vawsext,as.double(y),
                     as.logical(mask),
                     as.integer(nv),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     as.double(yext),
                     as.integer(nve),
                     hakt=as.double(hakt),
                     as.double(lambda0),
                     as.double(zobj$theta),
                     as.double(si2),
                     bi=as.double(zobj$bi),
                     theta=double(nv*n),
                     thext=double(nve*n),
                     as.integer(mc.cores),
                     double(prod(dlw)),
                     as.double(wghts),
                     double(nv*mc.cores),
                     double(nve*mc.cores))[c("bi","theta","thext","hakt")]
    dim(zobj$theta)<-c(nv,dy)
    thest <- array(zobj$thext,dyext)
    if(maxni) bi <- zobj$bi <- pmax(bi,zobj$bi)
    dim(zobj$bi)<-dy
    if(!is.null(u)) {
      cat("bandwidth: ",signif(hakt,3),"   MSE: ",
          signif(mean((zobj$theta-u)^2),3),"   MAE: ",
          signif(mean(abs(zobj$theta-u)),3)," mean(bi)=",
          signif(mean(zobj$bi),3),"\n")
      mae<-c(mae,signif(mean(abs(zobj$theta-u)),3))
    }
    if (max(total) >0) {
      cat(signif(total[k],2)*100,"% . ",sep="")
    }
  } else {
    thest <- NULL
  }
  cat("\n")
  list(y=y,theta=zobj$theta,yextsmoothed=thest,hakt=hakt,lambda=lambda,
       alpha=alpha,args=args,wghts=wghts,mae=mae,ni=zobj$bi)
}



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
                      data=NULL) {
  ##
  ##  this is the version with full size invcov (triangular storage)
  ##  needed for MPM
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
    cat("Corresponding bandwiths for specified correlation:",
        h0,
        "\n")
  }
  dim(invcov) <- c(nvec*nvec,n)
  hseq <- 1
  zobj <- list(bi = rep(1, n), theta = y)
  bi <- zobj$bi
  cat("Progress:")
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
    cat("step", k, "hakt", hakt, "time", format(Sys.time()), "\n")
    hseq <- c(hseq, hakt)
    dlw <- (2 * trunc(hakt / c(1, wghts)) + 1)
    if(k==kstar & !is.null(data)){
      zobj <- .Fortran(C_pvawse,
        as.double(y),
        as.double(data), ## data to smooth additionally
        as.logical(mask),
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
        bi = double(n), #binn
        theta = double(nvec * n),
        data = double(nsample*n),
        as.double(invcov[indcov,]),#compact storage
        as.integer(mc.cores),
        as.double(spmin),
        double(prod(dlw)),
        as.double(wghts),
        double(nvec * mc.cores),
        double(nsample * mc.cores),
        as.integer(np1),
        as.integer(np2),
        as.integer(np3))[c("bi", "theta", "hakt","data")]
        dim(zobj$data) <- dim(data)
    } else {
    zobj <- .Fortran(C_pvaws2,
      as.double(y),
      as.logical(mask),
      as.integer(nvec),
      as.integer(nvec * (nvec + 1) / 2),
      as.integer(n1),
      as.integer(n2),
      as.integer(n3),
      hakt = as.double(hakt),
      as.double(lambda0),
      as.double(zobj$theta),
      as.double(zobj$bi),
      bi = double(n), #binn
      theta = double(nvec * n),
      as.double(invcov[indcov,]),#compact storage
      as.integer(mc.cores),
      as.double(spmin),
      double(prod(dlw)),
      as.double(wghts),
      double(nvec * mc.cores),
      as.integer(np1),
      as.integer(np2),
      as.integer(np3))[c("bi", "theta", "hakt")]
    }
    dim(zobj$theta) <- c(nvec, dy)
    if (maxni)
      bi <- zobj$bi <- pmax(bi, zobj$bi)
    dim(zobj$bi) <- dy
    x <- 1.25 ^ k
    scorrfactor <- x / (3 ^ d * prod(scorr) * prod(h0) + x)
    lambda0 <- lambda * scorrfactor
    if (max(total) > 0) {
      cat(signif(total[k], 2) * 100, "% . ", sep = "")
    }
    k <- k + 1
    gc()
  }
  cat("\n")
  list(
#    y=y,
    theta=zobj$theta,
#    invcov=sweep(invcov, 2:(d + 1), zobj$bi, "*"),
    hakt=hakt,
    lambda=lambda,
    hseq = hseq,
    ni = zobj$bi,
    data= if(is.null(data)) data else zobj@data
  )
}

vpawscovm <- function(y,
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
                      data=NULL) {
  ##
  ##  this is the version with full size invcov (triangular storage)
  ##  needed for MPM
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
    cat("Corresponding bandwiths for specified correlation:",
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
  cat("Progress:")
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
    cat("step", k, "hakt", hakt, "time", format(Sys.time()), "\n")
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
        data[,mask] <- zobj@data
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
    if (max(total) > 0) {
      cat(signif(total[k], 2) * 100, "% . ", sep = "")
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
  cat("\n")
  list(
#    y=y,
    theta=theta,
#    invcov=sweep(invcov, 3:(d + 2), bi, "*"),
    hakt=hakt,
    lambda=lambda,
    hseq = hseq,
    ni = bi,
    data= if(is.null(data)) data else zobj@data
  )
}
