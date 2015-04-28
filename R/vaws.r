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
    zobj <- .Fortran("vaws",as.double(y),
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
                     double(nv*mc.cores),
                     PACKAGE="qMRI")[c("bi","theta","hakt")]
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
    zobj <- .Fortran("vawsext",as.double(y),
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
                     double(nve*mc.cores),
                     PACKAGE="qMRI")[c("bi","theta","thext","hakt")]
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

vsegm <- function(theta,bi,mask,level){
  dth <- dim(theta)
  dy <- dth[-1]
  nv <- dth[1]
  n1 <- dy[1]
  n2 <- dy[2]
  n3 <- dy[3]
  n <- n1*n2*n3
  segm <- .Fortran("vsegmen0",
                   as.double(theta),
                   as.double(bi),
                   as.logical(mask),
                   as.double(level),
                   double(prod(2*dy+1)),
                   as.integer(nv),
                   as.integer(2*n1+1),
                   as.integer(2*n2+1),
                   as.integer(2*n3+1),
                   integer(n),
                   segm=integer(n),
                   PACKAGE="aws")$segm
  dim(segm) <- dy
  segm          
}

fillsegm <- function(segm,mask,slevel){
  ##
  ##  replace segment number>slevel with a segment number 
  ##  from an adjacent region (in mask and segm<=slevel)
  ##
  ds <- dim(segm)
  n1 <- ds[1]
  n2 <- ds[2]
  n3 <- ds[3]
  n <- n1*n2*n3
  ntbl <- sum(segm[mask]>slevel)
  z <- .Fortran("fillsegm",
                nseg=as.integer(segm),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                as.logical(mask),
                integer(3*ntbl),
                as.integer(ntbl),
                as.integer(slevel),
                PACKAGE="aws")$nseg
  array(z,ds)
}

vdetrend <- function(theta,mask,h){
  ##
  ##  remove spatial trends using a median filter
  ##  trend removal will be everywhere as long as there are
  ##  active (mask==TRUE) voxel within a ball of radius h 
  ##
  dth <- dim(theta)
  dy <- dth[-1]
  nv <- dth[1]
  n1 <- dy[1]
  n2 <- dy[2]
  n3 <- dy[3]
  n <- n1*n2*n3
  newtheta <- theta
  mc.cores <- setCores(,reprt=FALSE)
  nwmd <- (2*as.integer(h)+1)^3
  parammd <- .Fortran("paramw3",
                      as.double(h),
                      as.double(c(1,1)),
                      ind=integer(3*nwmd),
                      w=double(nwmd),
                      n=as.integer(nwmd),
                      PACKAGE = "aws")[c("ind","n")]
  nwmd <- parammd$n
  parammd$ind <- parammd$ind[1:(3*nwmd)]
  dim(parammd$ind) <- c(3,nwmd)
  nind <- (2*h+1)^3
  for(k in 1:nv){
    z <- .Fortran("medsm1",
                  as.double(theta[k,,,]),
                  as.logical(mask),
                  as.integer(n1),
                  as.integer(n2),
                  as.integer(n3),
                  as.integer(parammd$ind),
                  as.integer(nwmd),
                  double(2*nwmd*mc.cores),
                  as.integer(mc.cores),
                  thnew=double(n),
                  PACKAGE="aws")$thnew
    newtheta[k,,,] <- z
  }
  theta-newtheta
}
gethani <- function(x,y,lkern,value,wght,eps=1e-2){
  .Fortran("gethani",
           as.double(x),
           as.double(y),
           as.integer(lkern),
           as.double(value),
           as.double(wght),
           as.double(eps),
           bw=double(1),
           PACKAGE="qMRI")$bw
}
