estimateIRfluid <- function(IRdata, InvTimes, segments,
                            TEScale = 100,
                            dataScale = 1000,
                            method = c("NLR", "QL"),
                            sigma = NULL,
                            L = 1,
                            varest = c("RSS","data"),
                            verbose = TRUE,
                            lower=c(0,0),
                            upper=c(2,2)){
   mask <- segments==1
   nvoxel <- sum(mask)
   ntimes <- length(InvTimes)
   itmax <- order(InvTimes)[ntimes]
   InvTimes[InvTimes==Inf] <- 50*max(InvTimes[InvTimes!=Inf])
   dimdata <- dim(IRdata)
   if(dimdata[1]!=ntimes) stop("estimateIRfluid: incompatible length of InvTimes")
   if(any(dimdata[-1]!=dim(mask))) stop("estimateIRfluid: incompatible dimension of segments")
   InvTimesScaled <- InvTimes/TEScale
   ## create necessary arrays
   npar <- 2 #  th2 for R, th1 for S
   Rx <- Sx <- Conv <- array(0,dim(mask))
   isConv <- array(0, nvoxel)
   isThresh <- array(FALSE, nvoxel)
   modelCoeff <- array(0, c(npar, nvoxel))
   if(varest[1]=="data"){
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
   if (method[1] == "QL") {
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
     order1 <- function(x) order(x)[1]
     itmin <- apply(IRdataFluid,2,order1)
     thetas[1,] <- IRdataFluid[itmax,]/dataScale
     thetas[2,] <- log(2)/InvTimes[itmin]*TEScale
     if (verbose){
        cat("Start estimation in", nvoxel, "voxel at", format(Sys.time()), "\n")
        pb <- txtProgressBar(0, nvoxel, style = 3)
     }
     for(xyz in 1:nvoxel){
     
     ivec <- IRdataFluid[, xyz]/dataScale
     th <- thetas[, xyz]
     
       res <- if (method[1] == "NLR") try(nls(ivec ~ IRhomogen(par, InvTimesScaled),
                                           data = list(InvTimesScaled),
                                           start = list(par = th),
                                           control = list(maxiter = 200,
                                                          warnOnly = TRUE)),silent=TRUE)
       else try(nls(ivec ~ IRhomogenQL(par, InvTimesScaled, CL, sig, L),
                    data = list(InvTimesScaled,
                                CL = CL,
                                sig = sig,
                                L = L),
                    start = list(par = th),
                    control = list(maxiter = 200,
                                   warnOnly = TRUE)),silent=TRUE)
       if (class(res) != "try-error"){
          thhat <- coef(res)
          outofrange <- any(thhat != pmin(upper,pmax(lower,thhat)))
       }
       if (class(res) == "try-error" || outofrange){
         # retry with port algorithm and bounds
         th <- pmin(upper,pmax(lower,th))
         res <- if (method[1] == "NLR") try(nls(ivec ~ IRhomogen(par, InvTimesScaled),
                                             data = list(InvTimes=InvTimesScaled),
                                             start = list(par = th),
                                               algorithm="port",
                                             control = list(maxiter = 200,
                                                            warnOnly = TRUE),
                                             lower=lower, upper=upper),silent=TRUE)
         else try(nls(ivec ~ IRhomogenQL(par, InvTimesScaled, CL, sig, L),
                      data = list(InvTimesScaled=InvTimesScaled,
                                  CL = CL,
                                  sig = sig,
                                  L = L),
                      start = list(par = th),
                      algorithm="port",
                      control = list(maxiter = 200,
                                     warnOnly = TRUE),
                      lower=lower, upper=upper),silent=TRUE)
       }
       if (class(res) != "try-error") {
         sres <- if(varest[1]=="RSS") getnlspars(res) else
           getnlspars2(res, shat[, xyz], sind )
         isConv[xyz] <- as.integer(res$convInfo$isConv)
         modelCoeff[, xyz] <- sres$coefficients
       }
       if (verbose) if(xyz%/%1000*1000==xyz) setTxtProgressBar(pb, xyz)
     }
  Rx[mask] <- modelCoeff[2,]
  Sx[mask] <- modelCoeff[1,]
  Conv[mask] <- isConv
  Sf <- median(modelCoeff[1,],na.rm=TRUE)
  Rf <- median(modelCoeff[2,],na.rm=TRUE)
  if (verbose){
    close(pb)
    cat("Finished estimation", format(Sys.time()), "\n","Sf",Sf,"Rf",Rf,"\n")
  }
  # Results are currently scaled by TEscale (R) and Datascale (S)
  list(Sf=Sf,Rf=Rf,Sx=Sx,Rx=Rx,sigma=sigma,Conv=Conv)
}
   


estimateIRsolid <- function(IRdata, InvTimes, segments, Sfluid, Rfluid,
                            TEScale = 100,
                            dataScale = 1000,
                            method = c("NLR", "QL"),
                            sigma = NULL,
                            L = 1,
                            varest = c("RSS","data"),
                            verbose = TRUE,
                            lower=c(0,0,0),
                            upper=c(.95,2,2)){
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
   Convx[segments==1] <- 1
   # set ICovx for fluid as (numerically) diag(rep(Inf),3)
   ICovx[1,1,segments==1] <- 1e20
   ICovx[2,2,segments==1] <- 1e20
   ICovx[3,3,segments==1] <- 1e20
   isConv <- array(FALSE, nvoxel)
   isThresh <- array(FALSE, nvoxel)
   modelCoeff <- array(0, c(npar, nvoxel))
   invCov <- array(0, c(npar, npar, nvoxel))
   rsigma <- array(0, nvoxel)
   if (method[1] == "QL") {
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
      thetas[2,] <- pmin(upper[2],pmax(lower[2],10/median(InvTimesScaled)))
      thetas[1,] <- 0.3
      if (verbose){
         cat("Start estimation in", nvoxel, "voxel at", format(Sys.time()), "\n")
         pb <- txtProgressBar(0, nvoxel, style = 3)
      }
      th1 <- (1:8)/10
      th2 <- Rfluid*c(.5,.6,.7,.8,.9,1.1,1.2)
      th3 <- Sfluid*(1:9)/10
      for(xyz in 1:nvoxel){
         
         ivec <- IRdataSolid[, xyz]/dataScale
         th <- thetas[, xyz]
##
##   initialize using grid search and optim
##
         best <- LSIRmix2(th,ivec,InvTimesScaled,Sfluid,Rfluid)
         for(i in 1:8) for(j in 1:7) for(k in 1:9){
            z <- LSIRmix2(c(th1[i],th2[j],th3[k]),ivec,InvTimesScaled,Sfluid,Rfluid)
            if(z < best){
               best <- z 
               th <- c(th1[i],th2[j],th3[k])
            }
         }                
         th <- pmin(upper,pmax(lower,th))
                           res <- if (method[1] == "NLR") try(optim(th, LSIRmix2, LSIRmix2grad, 
                                                    Y=ivec, InvTimes=InvTimesScaled, S0f=Sfluid, Rf=Rfluid,
                                                    method="L-BFGS-B",lower=lower,upper=upper))
                           else try(optim(th, LSIRmix2QL, LSIRmix2QLgrad, 
                                          Y=ivec, InvTimes=InvTimesScaled, S0f=Sfluid, Rf=Rfluid, 
                                          CL = CL, sig = sig, L = L,
                                          method="L-BFGS-B",lower=lower,upper=upper))
                  
                  if (class(res) != "try-error"){
           modelCoeff[,xyz] <- th <- res$par
           rsigma[xyz] <- sqrt(res$value)
           isConv[xyz] <- -res$convergence
         }
         res <- if (method[1] == "NLR") try(nls(ivec ~ IRmix2(par, ITS, Sfluid, Rfluid),
                                             data = list(ITS=InvTimesScaled, Sfluid=Sfluid, Rfluid=Rfluid),
                                             start = list(par = th),
                                             control = list(maxiter = 500,
                                                            warnOnly = TRUE)),silent=TRUE)
         else try(nls(ivec ~ IRmix2QL(par, ITS, Sfluid, Rfluid, CL, sig, L),
                      data = list(ITS=InvTimesScaled, Sfluid=Sfluid, Rfluid=Rfluid,
                                  CL = CL, sig = sig, L = L),
                      start = list(par = th),
                      control = list(maxiter = 500,
                                     warnOnly = TRUE)),silent=TRUE)
         if (class(res) != "try-error"){
           thhat <- coef(res)
           outofrange <- any(thhat != pmin(upper,pmax(lower,thhat)))
         }
         if (class(res) == "try-error" || outofrange){
            # retry with port algorithm and bounds
            th <- pmin(upper,pmax(lower,th))
            res <- if (method[1] == "NLR") try(nls(ivec ~ IRmix2(par, ITS, Sfluid, Rfluid),
                                                data = list(ITS=InvTimesScaled, Sfluid=Sfluid, Rfluid=Rfluid),
                                                start = list(par = th),
                                                algorithm="port",
                                                control = list(maxiter = 500,
                                                               warnOnly = TRUE),
                                                lower=lower, upper=upper),silent=TRUE)
            else try(nls(ivec ~ IRmix2QL(par, ITS, Sfluid, Rfluid, CL, sig, L),
                                     data = list(ITS=InvTimesScaled, Sfluid=Sfluid, Rfluid=Rfluid,
                                     CL = CL, sig = sig, L = L),
                         start = list(par = th),
                         algorithm="port",
                         control = list(maxiter = 500,
                                        warnOnly = TRUE),
                         lower=lower, upper=upper),silent=TRUE)
         }
         if (class(res) != "try-error") {
            sres <- if(varest[1]=="RSS") getnlspars(res) else
               getnlspars2(res, shat[, xyz], sind )
            isConv[xyz] <- as.integer(res$convInfo$stopCode)
            modelCoeff[, xyz] <- sres$coefficients
            if (sres$sigma != 0) {
               invCov[, , xyz] <- sres$invCov
               rsigma[xyz] <- sres$sigma
            }
         }
         if (verbose) if(xyz%/%1000*1000==xyz) setTxtProgressBar(pb, xyz)
      }
      if (verbose){
        close(pb)
        cat("Finished estimation", format(Sys.time()), "\n")
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


estimateIRsolid2 <- function(IRdata, InvTimes, segments, Sfluid, Rfluid,
                            TEScale = 100,
                            dataScale = 1000,
                            method = c("NLR", "QL"),
                            sigma = NULL,
                            L = 1,
                            varest = c("RSS","data"),
                            verbose = TRUE,
                            lower=c(0,0,0),
                            upper=c(.95,2,2)){
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
   df <- length(InvTimes)-npar
   fx <- Rx <- Sx <- rsdx <- array(0,dim(mask))
   ICovx <- array(0,c(3,3,prod(dim(mask))))
   Convx <- array(0,dim(mask))
   fx[segments==1] <- 1
   Rx[segments==1] <- Rfluid
   Sx[segments==1] <- Sfluid
   Convx[segments==1] <- 1
   # set ICovx for fluid as (numerically) diag(rep(Inf),3)
   ICovx[1,1,segments==1] <- 1e20
   ICovx[2,2,segments==1] <- 1e20
   ICovx[3,3,segments==1] <- 1e20
   isConv <- array(FALSE, nvoxel)
   isThresh <- array(FALSE, nvoxel)
   modelCoeff <- array(0, c(npar, nvoxel))
   invCov <- array(0, c(npar, npar, nvoxel))
   rsigma <- array(0, nvoxel)
   if (method[1] == "QL") {
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
   thetas[2,] <- pmin(upper[2],pmax(lower[2],10/median(InvTimesScaled)))
   thetas[1,] <- 0.3
   if (verbose){
      cat("Start estimation in", nvoxel, "voxel at", format(Sys.time()), "\n")
      pb <- txtProgressBar(0, nvoxel, style = 3)
   }
   th1 <- (1:8)/10
   th2 <- Rfluid*c(.5,.6,.7,.8,.9,1.1,1.2)
   th3 <- Sfluid*(1:9)/10
   for(xyz in 1:nvoxel){
      
      ivec <- IRdataSolid[, xyz]/dataScale
      th <- thetas[, xyz]
      ##
      ##   initialize using grid search and optim
      ##
      best <- LSIRmix2(th,ivec,InvTimesScaled,Sfluid,Rfluid)
      for(i in 1:8) for(j in 1:7) for(k in 1:9){
         z <- LSIRmix2(c(th1[i],th2[j],th3[k]),ivec,InvTimesScaled,Sfluid,Rfluid)
         if(z < best){
            best <- z 
            th <- c(th1[i],th2[j],th3[k])
         }
      }                
      th <- pmin(upper,pmax(lower,th))
      res <- if (method[1] == "NLR") try(optim(th, LSIRmix2, LSIRmix2grad, 
                                               Y=ivec, InvTimes=InvTimesScaled, S0f=Sfluid, Rf=Rfluid,
                                               method="L-BFGS-B",lower=lower,upper=upper,hessian=TRUE))
      else try(optim(th, LSIRmix2QL, LSIRmix2QLgrad, 
                     Y=ivec, InvTimes=InvTimesScaled, S0f=Sfluid, Rf=Rfluid, 
                     CL = CL, sig = sig, L = L,
                     method="L-BFGS-B",lower=lower,upper=upper,hessian=TRUE))
      
      if (class(res) != "try-error"){
         modelCoeff[,xyz] <- th <- res$par
         rsigma[xyz] <- sqrt(res$value/df)
         isConv[xyz] <- res$convergence
         invCov[,,xyz] <- res$hessian/res$value*df
      }
      if (verbose) if(xyz%/%1000*1000==xyz) setTxtProgressBar(pb, xyz)
   }
   if (verbose){
      close(pb)
      cat("Finished estimation", format(Sys.time()), "\n")
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
                                 varest = c("RSS","data"),
                                 verbose = TRUE,
                                 lower=c(0.05),
                                 upper=c(0.95)){
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
   if (method[1] == "QL") {
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
   thetas <- rep(0.3,nvoxel)
   if (verbose){
      cat("Start estimation in", nvoxel, "voxel at", format(Sys.time()), "\n")
      pb <- txtProgressBar(0, nvoxel, style = 3)
   }
   for(xyz in 1:nvoxel){
      
      ivec <- IRdataSolid[, xyz]/dataScale
      th <- thetas[xyz]
      Rs <- Rsm[xyz]
      Ss <- Ssm[xyz]
      
      res <- if (method[1] == "NLR") try(nls(ivec ~ IRmix2fix(par, ITS, Sf, Ss, Rf, Rs),
                                          data = list(ITS=InvTimesScaled, Sf=Sfluid, Ss=Ss, Rf=Rfluid, Rs=Rs),
                                          start = list(par = th),
                                          control = list(maxiter = 200,
                                                         warnOnly = TRUE)),silent=TRUE)
      else try(nls(ivec ~ IRmix2fixQL(par, ITS, Sf, Ss, Rf, Rs, CL, sig, L),
                   data = list(ITS=InvTimesScaled, Sf=Sfluid, Ss=Ss, Rf=Rfluid, Rs=Rs,
                               CL = CL, sig = sig, L = L),
                   start = list(par = th),
                   control = list(maxiter = 200,
                                  warnOnly = TRUE)),silent=TRUE)
      if (class(res) == "try-error"){
         # retry with port algorithm and bounds
         th <- pmin(upper,pmax(lower,th))
         res <- if (method[1] == "NLR") try(nls(ivec ~ IRmix2fix(par, ITS, Sf, Ss, Rf, Rs),
                                             data = list(ITS=InvTimesScaled, Sf=Sfluid, Ss=Ss, Rf=Rfluid, Rs=Rs),                                             start = list(par = th),
                                             algorithm="port",
                                             control = list(maxiter = 200,
                                                            warnOnly = TRUE),
                                             lower=lower, upper=upper),silent=TRUE)
         else try(nls(ivec ~ IRmix2fixQL(par, ITS, Sf, Ss, Rf, Rs, CL, sig, L),
                      data = list(ITS=InvTimesScaled, Sf=Sfluid, Ss=Ss, Rf=Rfluid, Rs=Rs,
                                  CL = CL, sig = sig, L = L),
                      start = list(par = th),
                      algorithm="port",
                      control = list(maxiter = 200,
                                     warnOnly = TRUE),
                      lower=lower, upper=upper),silent=TRUE)
      }
      if (class(res) != "try-error") {
         sres <- if(varest[1]=="RSS") getnlspars(res) else
            getnlspars2(res, shat[, xyz], sind )
         isConv[xyz] <- as.integer(res$convInfo$isConv)
         modelCoeff[xyz] <- sres$coefficients
         if (sres$sigma != 0) {
            invCov[ xyz] <- sres$invCov
            rsigma[xyz] <- sres$sigma
         }
      }
      if (verbose) if(xyz%/%1000*1000==xyz) setTxtProgressBar(pb, xyz)
   }
   if (verbose){
      close(pb)
      cat("Finished estimation", format(Sys.time()), "\n")
   }
fx[mask] <- pmin(upper,pmax(lower,modelCoeff))
ICovx[mask] <- invCov
Convx[mask] <- isConv
rsdx[mask] <- rsigma
# Results are currently scaled by TEscale (R) and Datascale (S)
list(fx=fx,Rx=Rsolid,Sx=Ssolid,Sf=Sfluid,Rf=Rfluid,ICovx=ICovx,Convx=Convx,sigma=sigma,rsdx=rsdx)
}

smoothIRSolid <- function(ergs,segm,kstar=24,ladjust=1){
   mask <- segm>1
   nvoxel <- sum(mask)
   bpars <- array(0,c(2,nvoxel))
   icovbpars <- array(0,c(2,2,nvoxel))
   bpars[1,] <- ergs$Rx[mask]
   bpars[2,] <- ergs$Sx[mask]
   ICovx <- ergs$ICovx
   dim(ICovx) <- c(3,3,prod(dim(mask)))
   icovbpars <- ICovx[-1,-1,mask]
   z <- vpawscov2(bpars, kstar, icovbpars/ladjust, segm>1)
   ergs$Rx[mask] <- z$theta[1,]
   ergs$Sx[mask] <- z$theta[2,]
   bi <- array(0,dim(mask))
   bi[mask] <- z$bi 
   ergs$bi <- bi
   ergs
   
}

estimateIR <- function(IRdata, InvTimes, segments, fixed=TRUE, smoothMethod=c("PAWS","Depth"),bw=5,
                       TEScale = 100,
                       dataScale = 1000,
                       method = c("NLR", "QL"),
                       sigma = NULL,
                       L = 1,
                       varest = c("RSS","data"),
                       kstar = 24,
                       ladjust = 1,
                       verbose = TRUE){
  
   ergsFluid <- estimateIRfluid(IRdata, InvTimes, segments)
   Sfluid <- median(ergsFluid$Sfluid)
   Rfluid <- median(ergsFluid$Rfluid)
   ergsBrain <- erstimateIRsolid(IRdata, InvTimes, segments, Sfluid, Rfluid)
   if(fixed) {
      if(smmothMethod[1]=="Depth") ergsSmooth <- SdepthSmooth(ergsBrain, segments)
      if(smmothMethod[1]=="PAWS") ergsSmooth <- smoothIRSolid(ergsBrain, segments, kstar, ladjust)
   }
   ergsSmooth
}

