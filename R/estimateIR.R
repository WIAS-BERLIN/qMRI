estimateIRfluid <- function(IRdata, InvTimes, segments,
                            TEScale = 100,
                            dataScale = 1000,
                            method = c("NLR", "QL"),
                            sigma = NULL,
                            L = 1,
                            maxR2star=50,
                            varest = c("RSS","data"),
                            verbose = TRUE,
                            lower=c(0,0),
                            upper=c(1,1)){
   mask <- segments==1
   nvoxel <- sum(mask)
   ntimes <- length(InvTimes)
   InvTimes[InvTimes==Inf] <- 50*max(InvTimes[InvTimes!=Inf])
   dimdata <- dim(IRdata)
   if(dimdata[1]!=ntimes) stop("estimateIRfluid: incompatible length of InvTimes")
   if(any(dimdata[-1]!=dim(mask))) stop("estimateIRfluid: incompatible dimension of segments")
   InvTimesScaled <- InvTimes/TEScale
   ## create necessary arrays
   npar <- 2 #  th2 for R, th1 for S
   Rx <- Sx <- array(0,dim(mask))
   isConv <- array(FALSE, nvoxel)
   isThresh <- array(FALSE, nvoxel)
   modelCoeff <- array(0, c(npar, nvoxel))
   if(varest=="data"){
     if(verbose) cat("estimating variance maps from data\n")
     ind <- (InvTimes == max(InvTimes))[1]
     ddata <- IRdata[ind,,,]
     shat <- aws::awsLocalSigma(ddata, steps=16,
                                        mask=(segments==1), ncoils=1, hsig=2.5,
                                        lambda=6,family="Gauss")$sigma
     dim(shat) <- dimdata[-1]
     shat <- shat[segments==1]
     shat[shat==0] <- quantile(shat,.8)
     shat <- shat
     if(is.null(sigma)) sigma <- median(shat) else shat <- NULL
   }
   if (method == "QL") {
   if(is.null(sigma)){ 
      method <- "NLR"
      warning("estimateIRfluid: method QL needs sigma estimated or supplied")
   }
     sig <- sigma/dataScale
     CL <- sig * sqrt(pi/2) * gamma(L + 0.5)/gamma(L)/gamma(1.5)
   }
   # initial parameters
     dim(IRdata) <- c(dimdata[1],prod(dim(segments)))
     IRdataFluid <- IRdata[,segments==1]
     thetas <- matrix(0,2,nvoxel)
     thetas[1,] <- IRdataFluid[(1:ntimes)[InvTimes == max(InvTimes)][1],]/dataScale
     thetas[2,] <- 1/median(InvTimesScaled)
     if (verbose){
        cat("Start estimation in", nvoxel, "voxel at", format(Sys.time()), "\n")
        pb <- txtProgressBar(0, nvoxel, style = 3)
     }
     for(xyz in 1:nvoxel){
     
     ivec <- IRdataFluid[, xyz]/dataScale
     th <- thetas[, xyz]
     
       res <- if (method == "NLR") try(nls(ivec ~ IRhomogen(par, InvTimesScaled),
                                           data = list(InvTimesScaled),
                                           start = list(par = th),
                                           control = list(maxiter = 200,
                                                          warnOnly = TRUE)))
       else try(nls(ivec ~ IRhomogenQL(par, InvTimesScaled, CL, sig, L),
                    data = list(InvTimesScaled,
                                CL = CL,
                                sig = sig,
                                L = L),
                    start = list(par = th),
                    control = list(maxiter = 200,
                                   warnOnly = TRUE)))
       if (class(res) == "try-error"){
         # retry with port algorithm and bounds
         th <- pmin(upper,pmax(lower,th))
         res <- if (method == "NLR") try(nls(ivec ~ IRhomogen(par, InvTimesScaled),
                                             data = list(InvTimes=InvTimesScaled),
                                             start = list(par = th),
                                               algorithm="port",
                                             control = list(maxiter = 200,
                                                            warnOnly = TRUE),
                                             lower=lower, upper=upper))
         else try(nls(ivec ~ IRhomogenQL(par, InvTimesScaled, CL, sig, L),
                      data = list(InvTimesScaled=InvTimesScaled,
                                  CL = CL,
                                  sig = sig,
                                  L = L),
                      start = list(par = th),
                      algorithm="port",
                      control = list(maxiter = 200,
                                     warnOnly = TRUE),
                      lower=lower, upper=upper))
       }
       if (class(res) != "try-error") {
         sres <- if(varest=="RSS") getnlspars(res) else
           getnlspars2(res, shat[, xyz], sind )
         isConv[xyz] <- as.integer(res$convInfo$isConv)
         modelCoeff[, xyz] <- sres$coefficients
       }
     }
  Rx[mask] <- modelCoeff[2,]
  Sx[mask] <- modelCoeff[1,]
  
  Sf <- median(modelCoeff[1,isConv==0])
  Rf <- median(modelCoeff[2,isConv==0])
  # Results are currently scaled by TEscale (R) and Datascale (S)
  list(Sf=Sf,Rf=Rf,Sx=Sx,Rx=Rx,sigma=sigma)
}
   


estimateIRsolid <- function(IRdata, InvTimes, segments, Sfluid, Rfluid,
                            TEScale = 100,
                            dataScale = 1000,
                            method = c("NLR", "QL"),
                            sigma = NULL,
                            L = 1,
                            maxR2star=50,
                            varest = c("RSS","data"),
                            verbose = TRUE,
                            lower=c(0,0),
                            upper=c(1,1)){
   mask <- segments>1
   nvoxel <- sum(mask)
   ntimes <- length(InvTimes)
   InvTimes[InvTimes==Inf] <- 50*max(InvTimes[InvTimes!=Inf])
   dimdata <- dim(IRdata)
   if(dimdata[1]!=ntimes) stop("estimateIRsolid: incompatible length of InvTimes")
   if(any(dimdata[-1]!=dim(mask))) stop("estimateIRsolid: incompatible dimension of segments")
   InvTimesScaled <- InvTimes/TEScale
   ## create necessary arrays
   npar <- 3 # th1 for f, th2 for R, th3 for S
   fx <- Rx <- Sx <- rsdx <- array(0,dim(mask))
   ICovx <- array(0,c(3,3,prod(dim(mask))))
   Convx <- array(0,dim(mask))
   fx[segments==1] <- 1
   Rx[segments==1] <- Rfluid
   Sx[segments==1] <- Sfluid
   Convx[segments==1] <- 0
   # set ICovx for fluid as (numerically) diag(rep(Inf),3)
   ICovx[1,1,segments==1] <- 1e20
   ICovx[2,2,segments==1] <- 1e20
   ICovx[3,3,segments==1] <- 1e20
   isConv <- array(FALSE, nvoxel)
   isThresh <- array(FALSE, nvoxel)
   modelCoeff <- array(0, c(npar, nvoxel))
   invCov <- array(0, c(npar, npar, nvoxel))
   rsigma <- array(0, nvoxel)
   if (method == "QL") {
      if(is.null(sigma)){ 
         method <- "NLR"
         warning("estimateIRsolid: method QL needs sigma estimated from fluid or supplied")
      }
         sig <- sigma/dataScale
         CL <- sig * sqrt(pi/2) * gamma(L + 0.5)/gamma(L)/gamma(1.5)
      }
      # initial parameters
      dim(IRdata) <- c(dimdata[1],prod(dim(segments)))
      IRdataSolid <- IRdata[,mask]
      thetas <- matrix(0,3,nvoxel)
      thetas[3,] <- IRdataSolid[(1:ntimes)[InvTimes == max(InvTimes)][1],]/dataScale
      thetas[2,] <- 1/median(InvTimesScaled)
      thetas[1,] <- 0.3
      if (verbose){
         cat("Start estimation in", nvoxel, "voxel at", format(Sys.time()), "\n")
         pb <- txtProgressBar(0, nvoxel, style = 3)
      }
      for(xyz in 1:nvoxel){
         
         ivec <- IRdataSolid[, xyz]/dataScale
         th <- thetas[, xyz]
         
         res <- if (method == "NLR") try(nls(ivec ~ IRmix2(par, ITS, Sfluid, Rfluid),
                                             data = list(ITS=InvTimesScaled, Sfluid=Sfluid, Rfluid=Rfluid),
                                             start = list(par = th),
                                             control = list(maxiter = 200,
                                                            warnOnly = TRUE)))
         else try(nls(ivec ~ IRmix2QL(par, ITS, Sfluid, Rfluid, CL, sig, L),
                      data = list(ITS=InvTimesScaled, Sfluid=Sfluid, Rfluid=Rfluid,
                                  CL = CL, sig = sig, L = L),
                      start = list(par = th),
                      control = list(maxiter = 200,
                                     warnOnly = TRUE)))
         if (class(res) == "try-error"){
            # retry with port algorithm and bounds
            th <- pmin(upper,pmax(lower,th))
            res <- if (method == "NLR") try(nls(ivec ~ IRmix2(par, ITS, Sfluid, Rfluid),
                                                data = list(ITS=InvTimesScaled, Sfluid=Sfluid, Rfluid=Rfluid),
                                                start = list(par = th),
                                                algorithm="port",
                                                control = list(maxiter = 200,
                                                               warnOnly = TRUE),
                                                lower=lower, upper=upper))
            else try(nls(ivec ~ IRmix2QL(par, ITS, Sfluid, Rfluid, CL, sig, L),
                                     data = list(ITS=InvTimesScaled, Sfluid=Sfluid, Rfluid=Rfluid,
                                     CL = CL, sig = sig, L = L),
                         start = list(par = th),
                         algorithm="port",
                         control = list(maxiter = 200,
                                        warnOnly = TRUE),
                         lower=lower, upper=upper))
         }
         if (class(res) != "try-error") {
            sres <- if(varest=="RSS") getnlspars(res) else
               getnlspars2(res, shat[, xyz], sind )
            isConv[xyz] <- as.integer(res$convInfo$isConv)
            modelCoeff[, xyz] <- sres$coefficients
            if (sres$sigma != 0) {
               invCov[, , xyz] <- sres$invCov
               rsigma[xyz] <- sres$sigma
            }
         }
      }
      fx[mask] <- modelCoeff[1,]
      Rx[mask] <- modelCoeff[2,]
      Sx[mask] <- modelCoeff[3,]
      ICovx[,,mask] <- invCov
      dim(ICovx) <- c(3,3,dim(mask))
      Convx[mask] <- isConv
      rsdx[mask] <- rsigma
# Results are currently scaled by TEscale (R) and Datascale (S)
      list(fx=fx,Rx=Rx,Sx=Sx,Sf=Sfluid,Rf=Rfluid,ICovx=ICovx,Convx=Convx,sigma=sigma,rsdx=rsdx)
   }



estimateIRsolidfixed <- function(IRdata, InvTimes, segments, Sfluid, Rfluid, Ssolid, Rsolid,
                                 TEScale = 100,
                                 dataScale = 1000,
                                 method = c("NLR", "QL"),
                                 sigma = NULL,
                                 L = 1,
                                 maxR2star=50,
                                 varest = c("RSS","data"),
                                 verbose = TRUE,
                                 lower=c(0,0),
                                 upper=c(1,1)){
   mask <- segments>1
   nvoxel <- sum(mask)
   ntimes <- length(InvTimes)
   InvTimes[InvTimes==Inf] <- 50*max(InvTimes[InvTimes!=Inf])
   dimdata <- dim(IRdata)
   if(dimdata[1]!=ntimes) stop("estimateIRsolid: incompatible length of InvTimes")
   if(any(dimdata[-1]!=dim(mask))) stop("estimateIRsolid: incompatible dimension of segments")
   InvTimesScaled <- InvTimes/TEScale
   ## create necessary arrays
   npar <- 1 # th1 for f
   fx  <- rsdx <- array(0,dim(mask))
   ICovx <- array(0,prod(dim(mask)))
   Convx <- array(0,dim(mask))
   fx[segments==1] <- 1
   Rx <- Rsolid
   Sx <- Ssolid
   Convx[segments==1] <- 0
   # set ICovx for fluid as (numerically) diag(rep(Inf),3)
   ICovx[segments==1] <- 1e20
   isConv <- array(FALSE, nvoxel)
   isThresh <- array(FALSE, nvoxel)
   modelCoeff <- numeric(nvoxel)
   invCov <- numeric(nvoxel)
   rsigma <- numeric(nvoxel)
   if (method == "QL") {
      if(is.null(sigma)){ 
         method <- "NLR"
         warning("estimateIRsolid: method QL needs sigma estimated from fluid or supplied")
      }
      sig <- sigma/dataScale
      CL <- sig * sqrt(pi/2) * gamma(L + 0.5)/gamma(L)/gamma(1.5)
   }
   # initial parameters
   dim(IRdata) <- c(dimdata[1],prod(dim(segments)))
   IRdataSolid <- IRdata[,mask]
   Rsm <- Rsolid[mask]
   Ssm <- Ssolid[mask]
   thetas <- rep(0.1,nvoxel)
   if (verbose){
      cat("Start estimation in", nvoxel, "voxel at", format(Sys.time()), "\n")
      pb <- txtProgressBar(0, nvoxel, style = 3)
   }
   for(xyz in 1:nvoxel){
      
      ivec <- IRdataSolid[, xyz]/dataScale
      th <- thetas[, xyz]
      Rs <- Rsm[xyz]
      Ss <- Ssm[xyz]
      
      res <- if (method == "NLR") try(nls(ivec ~ IRmix2fix(par, ITS, Sf, Ss, Rf, Rs),
                                          data = list(ITS=InvTimesScaled, Sf=Sfluid, Ss=Ss, Rf=Rfluid, Rs=Rs),
                                          start = list(par = th),
                                          control = list(maxiter = 200,
                                                         warnOnly = TRUE)))
      else try(nls(ivec ~ IRmix2fixQL(par, ITS, Sf, Ss, Rf, Rs, CL, sig, L),
                   data = list(ITS=InvTimesScaled, Sf=Sfluid, Ss=Ss, Rf=Rfluid, Rs=Rs,
                               CL = CL, sig = sig, L = L),
                   start = list(par = th),
                   control = list(maxiter = 200,
                                  warnOnly = TRUE)))
      if (class(res) == "try-error"){
         # retry with port algorithm and bounds
         th <- pmin(upper,pmax(lower,th))
         res <- if (method == "NLR") try(nls(ivec ~ IRmix2fix(par, ITS, Sf, Ss, Rf, Rs),
                                             data = list(ITS=InvTimesScaled, Sf=Sfluid, Ss=Ss, Rf=Rfluid, Rs=Rs),                                             start = list(par = th),
                                             algorithm="port",
                                             control = list(maxiter = 200,
                                                            warnOnly = TRUE),
                                             lower=lower, upper=upper))
         else try(nls(ivec ~ IRmix2fixQL(par, ITS, Sf, Ss, Rf, Rs, CL, sig, L),
                      data = list(ITS=InvTimesScaled, Sf=Sfluid, Ss=Ss, Rf=Rfluid, Rs=Rs,
                                  CL = CL, sig = sig, L = L),
                      start = list(par = th),
                      algorithm="port",
                      control = list(maxiter = 200,
                                     warnOnly = TRUE),
                      lower=lower, upper=upper))
      }
      if (class(res) != "try-error") {
         sres <- if(varest=="RSS") getnlspars(res) else
            getnlspars2(res, shat[, xyz], sind )
         isConv[xyz] <- as.integer(res$convInfo$isConv)
         modelCoeff[xyz] <- sres$coefficients
         if (sres$sigma != 0) {
            invCov[, , xyz] <- sres$invCov
            rsigma[xyz] <- sres$sigma
         }
      }
   }
fx[mask] <- modelCoeff
ICovx[mask] <- invCov
Convx[mask] <- isConv
rsdx[mask] <- rsigma
# Results are currently scaled by TEscale (R) and Datascale (S)
list(fx=fx,Rx=Rx,Sx=Sx,Sf=Sfluid,Rf=Rfluid,ICovx=ICovx,Convx=Convx,sigma=sigma,rsdx=rsdx)
}

estimateIR <- function(IRdata, InvTimes, segments, fixed=TRUE, smoothMethod=c("Depth","PAWS"),bw=5,
                       TEScale = 100,
                       dataScale = 1000,
                       method = c("NLR", "QL"),
                       sigma = NULL,
                       L = 1,
                       maxR2star=50,
                       varest = c("RSS","data"),
                       verbose = TRUE){
  
   ergsFluid <- estimateIRfluid(IRdata, InvTimes, segments)
   Sfluid <- median(ergsFluid$Sfluid)
   Rfluid <- median(ergsFluid$Rfluid)
   ergsBrain <- erstimateIRsolid(IRdata, InvTimes, segments, Sfluid, Rfluid)
   if(fixed) {
      if(smmothMethod=="Depth") ergsSmooth <- SdepthSmooth(ergsBrain, segments)
      if(smmothMethod=="PAWS") ergsSmooth <- vpawsaws::vpawscov2(ergsBrain, segments)
   }
   ergsSmooth
}

