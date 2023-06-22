plmatrix <- function(x, FUN, ..., mc.cores = setCores(,reprt=FALSE)){
  dx <- dim(x)[2]
  if(mc.cores>dx) mc.cores <- dx
  n <- trunc((dx-1)/mc.cores)+1
  lx <- list(NULL)
  for(i in 1:(mc.cores-1)) lx[[i]] <- x[,(i-1)*n+1:n]
  lx[[mc.cores]] <- x[,((mc.cores-1)*n+1):dx]
  cl <- makeCluster(mc <- mc.cores)
  lz <- parLapply(cl, lx, FUN , ...)
  stopCluster(cl)
  z <- matrix(0,length(lz[[1]])/n, dx)
  for(i in 1:(mc.cores-1)) z[,(i-1)*n+1:n] <- lz[[i]]
  z[,((mc.cores-1)*n+1):dx] <- lz[[mc.cores]]
  z
}
# in: thetas , IRdataFluid, InvTimesScaled
# out: isConv, modelcoef
pIRfluid <- function(x, InvTimesScaled, method, sigma, CL, sig, L, varest, lower=c(0,0),
                     upper=c(2,2)){
   nvoxel <- dim(x)[2]
   npar <- 2
   ntimes <- length(InvTimesScaled)
   sind <- rep(1,ntimes)
   thetas <- x[1:2,]
   IRdataFluid <- x[-c(1:2),]
ergs <- array(0,c(npar+1,nvoxel))

for(xyz in 1:nvoxel){
  
  ivec <- IRdataFluid[, xyz]
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
  if (!inherits(res, "try-error")){
    thhat <- coef(res)
    outofrange <- any(thhat != pmin(upper,pmax(lower,thhat)))
  }
  if (inherits(res, "try-error") || outofrange){
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
  if (!inherits(res, "try-error")) {
    sres <- if(varest[1]=="RSS") getnlspars(res) else
      getnlspars2(res, sigma, sind )
    ergs[npar+1,xyz] <- as.integer(res$convInfo$isConv)
    ergs[1:npar, xyz] <- sres$coefficients
  }
}
ergs
}

pIRsolid <- function(x, InvTimesScaled, Rfluid, Sfluid, method, sigma, CL, sig, L,
                     varest, 
                     lower=c(0,0,0),
                     upper=c(.95,2,2)){
  nvoxel <- dim(x)[2]
  npar <- 3
  thetas <- x[1:npar,]
  IRdataSolid <- x[-c(1:npar),]
  ntimes <- length(InvTimesScaled)
  sind <- rep(1,ntimes)
  df <- ntimes-3
  ergs <- array(0,c(npar+npar*npar+2,nvoxel))
  th1 <- (1:8)/10
  th2 <- Rfluid*c(.5,.6,.7,.8,.9,1.1,1.2)
  th3 <- Sfluid*(1:9)/10
  for(xyz in 1:nvoxel){
    
    ivec <- IRdataSolid[, xyz]
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
    
    if (!inherits(res, "try-error")){
      ergs[1:npar,xyz] <- th <- res$par
      ergs[npar+npar*npar+1,xyz] <- sqrt(res$value)
      ergs[npar+npar*npar+2,xyz] <- -res$convergence
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
    if (!inherits(res, "try-error")){
      thhat <- coef(res)
      outofrange <- any(thhat != pmin(upper,pmax(lower,thhat)))
    }
    if (inherits(res, "try-error") || outofrange){
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
    if (!inherits(res, "try-error")) {
      sres <- if(varest[1]=="RSS") getnlspars(res) else
        getnlspars2(res, sigma, sind )
      ergs[npar+npar*npar+2,xyz] <- as.integer(res$convInfo$stopCode)
      ergs[1:npar, xyz] <- sres$coefficients
      if (sres$sigma != 0) {
        ergs[npar+1:(npar*npar) , xyz] <- as.vector(sres$invCov)*df
        ergs[npar+npar*npar+1,xyz] <- sres$sigma
      }
    }
  }  
  ergs
}
