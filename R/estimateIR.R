estimateIRfluid <- function(IRdata, InvTimes, segments,
                            TEScale = 100,
                            dataScale = 1000,
                            method = c("NLR", "QL"),
                            sigma = NULL,
                            L = 1,
                            maxR2star=50,
                            varest = c("RSS","data"),
                            verbose = TRUE){
   mask <- segments==1
   nvoxel <- sum(mask)
   ntimes <- length(InvTimes)
   InvTimes[InvTimes==Inf] <- 50*max(InvTimes[InvTimes!=Inf])
   dimdata <- dim(IRdata)
   if(dimdata[1]!=ntimes) stop("estimateIRfluid: incompatible length of InvTimes")
   if(any(dimdata[-1]!=dim(mask))) stop("estimateIRfluid: incompatible dimension of segments")
   InvTimesScaled <- InvTimes/TEScale
   ## create necessary arrays
   isConv <- array(FALSE, nvoxel)
   isThresh <- array(FALSE, nvoxel)
   modelCoeff <- array(0, c(npar, nvoxel))
   invCov <- array(0, c(npar, npar, nvoxel))
   rsigma <- array(0, nvoxel)
   if (verbose){
     cat("Start estimation in", nvoxel, "voxel at", format(Sys.time()), "\n")
     pb <- txtProgressBar(0, nvoxel, style = 3)
   }
   if(varest=="data"){
     if(verbose) cat("estimating variance maps from data\n")
     ind <- (InvTimes == max(InvTimes))[1]
     ddata <- IRdata[ind,,,]
     shat <- aws::awsLocalSigma(ddata, steps=16,
                                        mask=(segments==1), ncoils=1, hsig=2.5,
                                        lambda=6,family="Gauss")$sigma
     }
     dim(shat) <- dimdata[-1]
     shat <- shat[segments==1]
     shat[shat==0] <- quantile(shat,.8)
     shat <- shat
     if(is.null(sigma)) sigma <- median(shat)
   } else shat <- NULL
   if (method == "QL") {
     sig <- sigma/dataScale
     CL <- sigma * sqrt(pi/2) * gamma(L + 0.5)/gamma(L)/gamma(1.5)
   }
   # initial parameters
     IRdataFluid <- IRdata[,segments==1]
     th <- matrix(0,2,nvoxel)
     th[1,] <- IRdataFluid[(InvTimes == max(InvTimes))[1],]/dataScale
     th[2,] <- median(InvTimesScaled)
   for(xyz in 1:nvoxel){
     
     ivec <- IRdataFluid[, xyz]/dataScale
     th <- thetas[, xyz]
     
       res <- if (method == "NLR") try(nls(ivec ~ IRhomogen(par, InvTimes),
                                           data = list(InvTimesScaled),
                                           start = list(par = th),
                                           weights = wghts,
                                           control = list(maxiter = 200,
                                                          warnOnly = TRUE)))
       else try(nls(ivec ~ IRhomogenQL(par, InvTimes, CL, sig, L),
                    data = list(InvTimes,
                                CL = CL,
                                sig = sig,
                                L = L),
                    start = list(par = th),
                    weights = wghts,
                    control = list(maxiter = 200,
                                   warnOnly = TRUE)))
       if (class(res) == "try-error"){
         # retry with port algorithm and bounds
         th <- pmin(upper,pmax(lower,th))
         res <- if (method == "NLR") try(nls(ivec ~ IRhomogen(par, InvTimes),
                                             data = list(InvTimes=InvTimesScaled),
                                             start = list(par = th),
                                               algorithm="port",
                                             weights = wghts,
                                             control = list(maxiter = 200,
                                                            warnOnly = TRUE),
                                             lower=lower, upper=upper))
         else try(nls(ivec ~ IRhomogenQL(par, xmat, CL, sig, L),
                      data = list(InvTimes=InvTimesScaled,
                                  CL = CL,
                                  sig = sig,
                                  L = L),
                      start = list(par = th),
                      algorithm="port",
                      weights = wghts,
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
  Sf <- median(modelCoeff[1,isConv==0])
  Rf <- median(modelCoeff[2,isConv==0])
  list(Sf=Sf,Rf=Rf,sigma=sigma)
}
   


estimateIRsolid <- function(IRdata, InvTimes, segments, Sfluid, Tfluid,
                            TEScale = 100,
                            dataScale = 1000,
                            method = c("NLR", "QL"),
                            sigma = NULL,
                            L = 1,
                            maxR2star=50,
                            varest = c("RSS","data"),
                            verbose = TRUE){
  
}

estimateIRsolidfixed <- function(IRdata, InvTimes, segments, Sfluid, Tfluid, Ssolid, Tsolid,
                                 TEScale = 100,
                                 dataScale = 1000,
                                 method = c("NLR", "QL"),
                                 sigma = NULL,
                                 L = 1,
                                 maxR2star=50,
                                 varest = c("RSS","data"),
                                 verbose = TRUE){
  
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
   Tfluid <- median(ergsFluid$Tfluid)
   ergsBrain <- erstimateIRsolid(IRdata, InvTimes, segments, Sfluid, Tfluid)
   if(fixed) {
      if(smmothMethod=="Depth") ergsSmooth <- SdepthSmooth(ergsBrain, segments)
      if(smmothMethod=="PAWS") ergsSmooth <- vpawsaws::vpawscov2(ergsBrain, segments)
   }
}

