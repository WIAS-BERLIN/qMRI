estatics3 <- function(par, design){
  ##
  ## former: qflashpl
  ##
  ## ESTATICS model with 3+1 parameters
  ##
  ## S_{T1} = par[1] * exp(- par[4] * TE)
  ## S_{MT} = par[2] * exp(- par[4] * TE)
  ## S_{PD} = par[3] * exp(- par[4] * TE)
  ##
  n <- dim(design)[1]
  z <- .Fortran(C_estatics3,
                as.double(par),
                as.double(design),
                as.integer(n),
                fval = double(n),
                grad = double(4*n))[c("fval", "grad")]
  fval <- z$fval
  attr(fval,"gradient") <- matrix(z$grad, n, 4)
  fval
}

estatics2 <- function(par, design){
  ##
  ## former: qflashpl2
  ##
  ## ESTATICS model with 2+1 parameters
  ##
  ## S_{T1} = par[1] * exp(- par[3] * TE)
  ## S_{PD} = par[2] * exp(- par[3] * TE)
  ##
  n <- dim(design)[1]
  z <- .Fortran(C_estatics2,
                as.double(par),
                as.double(design),
                as.integer(n),
                fval = double(n),
                grad = double(3*n))[c("fval", "grad")]
  fval <- z$fval
  attr(fval, "gradient") <- matrix(z$grad, n, 3)
  fval
}

estatics1 <- function(par, design){
  ##
  ## former: qflashpl3
  ##
  ## "ESTATICS" model with 1+1 parameters
  ##
  ## S_{T1} = par[1] * exp(- par[2] * TE)
  ##
  n <- dim(design)[1]
  z <- .Fortran(C_estatics1,
                as.double(par),
                as.double(design),
                as.integer(n),
                fval = double(n),
                grad = double(3*n))[c("fval", "grad")]
  fval <- z$fval
  attr(fval, "gradient") <- matrix(z$grad, n, 3)
  fval
}

estatics3fixedR2 <- function(par, R2star, design){
  ##
  ## former: qflashpl0
  ##
  ## ESTATICS model with 3 parameters
  ##
  ## S_{T1} = par[1] * exp(- R2star * TE)
  ## S_{MT} = par[2] * exp(- R2star * TE)
  ## S_{PD} = par[3] * exp(- R2star * TE)
  ##
  n <- dim(design)[1]
  z <- .Fortran(C_estatics3fixedr2,
                as.double(par),
                as.double(R2star),
                as.double(design),
                as.integer(n),
                fval = double(n),
                grad = double(3*n))[c("fval", "grad")]
  fval <- z$fval
  attr(fval, "gradient") <- matrix(z$grad, n, 3)
  fval
}

estatics2fixedR2 <- function(par, R2star, design){
  ##
  ## former: qflashpl20
  ##
  ## ESTATICS model with 2 parameters
  ##
  ## S_{T1} = par[1] * exp(- R2star * TE)
  ## S_{PD} = par[2] * exp(- R2star * TE)
  ##
  n <- dim(design)[1]
  z <- .Fortran(C_estatics2fixedr2,
                as.double(par),
                as.double(R2star),
                as.double(design),
                as.integer(n),
                fval = double(n),
                grad = double(2*n))[c("fval", "grad")]
  fval <- z$fval
  attr(fval, "gradient") <- matrix(z$grad, n, 2)
  fval
}

estatics1fixedR2 <- function(par, R2star, design){
  ##
  ## ESTATICS model with 1 parameter
  ##
  ## S_{T1} = par[1] * exp(- R2star * TE)
  ##
  n <- dim(design)[1]
  z <- .Fortran(C_estatics1fixedr2,
                as.double(par),
                as.double(R2star),
                as.double(design),
                as.integer(n),
                fval = double(n),
                grad = double(2*n))[c("fval", "grad")]
  fval <- z$fval
  attr(fval, "gradient") <- matrix(z$grad, n, 2)
  fval
}

estatics3QL <- function(par, design, CL, sigma, L){
  ##
  ## former: qflashplQL
  ##
  ## ESTATICS model with 3+1 parameters with QL
  ##
  ## S_{T1} = par[1] * exp(- par[4] * TE)
  ## S_{MT} = par[2] * exp(- par[4] * TE)
  ## S_{PD} = par[3] * exp(- par[4] * TE)
  ##
  n <- dim(design)[1]
  z <- .Fortran(C_estatics3,
                as.double(par),
                as.double(design),
                as.integer(n),
                fval = double(n),
                grad = double(4*n))[c("fval", "grad")]
  sfval <- z$fval/sigma
  fval <- CL * hg1f1(-.5, L, -sfval*sfval/2)
  CC <- CL * hg1f1(.5, L+1, -sfval*sfval/2) * sfval /2/L/sigma
  grad <- matrix(CC*z$grad, n, 4)
  ind <- sfval>100
  if(any(ind)){
# use LS if there is no real difference, factor to keep monotonicity
     fval[ind] <- z$fval[ind]*1.0001
     grad[ind,] <- matrix(z$grad, n, 4)[ind, ]*1.0001
  }
  attr(fval, "gradient") <- grad
  fval
}

estatics2QL <- function(par, design, CL, sigma, L){
  ##
  ## former: qflashpl2QL
  ##
  ## ESTATICS model with 2+1 parameters with QL
  ##
  ## S_{T1} = par[1] * exp(- par[3] * TE)
  ## S_{PD} = par[2] * exp(- par[3] * TE)
  ##
  n <- dim(design)[1]
  z <- .Fortran(C_estatics2,
                as.double(par),
                as.double(design),
                as.integer(n),
                fval = double(n),
                grad = double(3*n))[c("fval", "grad")]
  # CL <- sigma * sqrt(pi/2) * gamma(L+0.5) / gamma(L) / gamma(1.5)
  sfval <- z$fval/sigma
  fval <- CL * hg1f1(-.5, L, -sfval*sfval/2)
  CC <- CL * hg1f1(.5, L+1, -sfval*sfval/2) * sfval /2/L/sigma
  grad <- matrix(CC*z$grad, n, 3)
  ind <- sfval>100
  if(any(ind)){
# use LS if there is no real difference, factor to keep monotonicity
    fval[ind] <- z$fval[ind]*1.0001
    grad[ind,] <- matrix(z$grad, n, 3)[ind, ]*1.0001
  }
  attr(fval, "gradient") <- grad
  fval
}

estatics1QL <- function(par, design, CL, sigma, L){
  ##
  ## "ESTATICS" model with 1+1 parameters with QL
  ##
  ## S_{T1} = par[1] * exp(- par[2] * TE)
  ##
  n <- dim(design)[1]
  z <- .Fortran(C_estatics1,
                as.double(par),
                as.double(design),
                as.integer(n),
                fval = double(n),
                grad = double(2*n))[c("fval", "grad")]
  # CL <- sigma * sqrt(pi/2) * gamma(L+0.5) / gamma(L) / gamma(1.5)
  sfval <- z$fval/sigma
  fval <- CL * hg1f1(-.5, L, -sfval*sfval/2)
  CC <- CL * hg1f1(.5, L+1, -sfval*sfval/2) * sfval /2/L/sigma
  grad <- matrix(CC*z$grad, n, 2)
  ind <- sfval>100
  if(any(ind)){
# use LS if there is no real difference, factor to keep monotonicity
    fval[ind] <- z$fval[ind]*1.0001
    grad[ind,] <- matrix(z$grad, n, 2)[ind, ]*1.0001
  }
  attr(fval, "gradient") <- grad
  fval
}

estatics3QLfixedR2 <- function(par, R2star, design, CL, sigma, L){
  ##
  ## former: qflashpl0QL
  ##
  ## ESTATICS model with 3 parameters
  ##
  ## S_{T1} = par[1] * exp(- R2star * TE)
  ## S_{MT} = par[2] * exp(- R2star * TE)
  ## S_{PD} = par[3] * exp(- R2star * TE)
  ##
  n <- dim(design)[1]
  z <- .Fortran(C_estatics3fixedr2,
                as.double(par),
                as.double(R2star),
                as.double(design),
                as.integer(n),
                fval = double(n),
                grad = double(3*n))[c("fval", "grad")]
  # CL <- sigma * sqrt(pi/2) * gamma(L+0.5) / gamma(L) / gamma(1.5)
  sfval <- z$fval/sigma
  fval <- CL * hg1f1(-.5, L, -sfval*sfval/2)
  CC <- CL * hg1f1(.5, L+1, -sfval*sfval/2) * sfval /2/L/sigma
  grad <- matrix(CC*z$grad, n, 3)
  ind <- sfval>100
  if(any(ind)){
# use LS if there is no real difference, factor to keep monotonicity
    fval[ind] <- z$fval[ind]*1.0001
    grad[ind,] <- matrix(z$grad, n, 3)[ind, ]*1.0001
  }
  attr(fval, "gradient") <- grad
  fval
}

estatics2QLfixedR2 <- function(par, R2star, design, CL, sigma, L){
  ##
  ##  former: qflashpl20QL
  ##
  ## ESTATICS model with 2 parameters
  ##
  ## S_{T1} = par[1] * exp(- R2star * TE)
  ## S_{MT} = par[2] * exp(- R2star * TE)
  ## S_{PD} = par[3] * exp(- R2star * TE)
  ##
  n <- dim(design)[1]
  z <- .Fortran(C_estatics2fixedr2,
                as.double(par),
                as.double(R2star),
                as.double(design),
                as.integer(n),
                fval = double(n),
                grad = double(2*n))[c("fval", "grad")]
  # CL <- sigma * sqrt(pi/2) * gamma(L+0.5) / gamma(L) / gamma(1.5)
  sfval <- pmin(z$fval/sigma,1e10)
  fval <- CL * hg1f1(-.5, L, -sfval*sfval/2)
  CC <- CL * hg1f1(.5, L+1, -sfval*sfval/2) * sfval /2/L/sigma
  grad <- matrix(CC*z$grad, n, 2)
  ind <- sfval>100
  if(any(ind)){
# use LS if there is no real difference, factor to keep monotonicity
    fval[ind] <- z$fval[ind]*1.0001
    grad[ind,] <- matrix(z$grad, n, 2)[ind, ]*1.0001
  }
  attr(fval, "gradient") <- grad
  fval
}

estatics1QLfixedR2 <- function(par, R2star, design, CL, sigma, L){
  ##
  ## ESTATICS model with 1 parameters
  ##
  ## S_{T1} = par[1] * exp(- R2star * TE)
  ##
  n <- dim(design)[1]
  z <- .Fortran(C_estatics1fixedr2,
                as.double(par),
                as.double(R2star),
                as.double(design),
                as.integer(n),
                fval = double(n),
                grad = double(n))[c("fval", "grad")]
  # CL <- sigma * sqrt(pi/2) * gamma(L+0.5) / gamma(L) / gamma(1.5)
  sfval <- pmin(z$fval/sigma,1e10)
  fval <- CL * hg1f1(-.5, L, -sfval*sfval/2)
  CC <- CL * hg1f1(.5, L+1, -sfval*sfval/2) * sfval /2/L/sigma
  grad <- matrix(CC*z$grad, n, 1)
  ind <- sfval>100
  if(any(ind)){
# use LS if there is no real difference, factor to keep monotonicity
    fval[ind] <- z$fval[ind]*1.0001
    grad[ind,] <- matrix(z$grad, n, 1)[ind, ]*1.0001
  }
  attr(fval, "gradient") <- grad
  fval
}

reparamtrizeconstTR <- function(th,a1,a2,TR){
  #
  #   Parameter transformation for constant TR
  #   from    S_T1(TE,alpha) = th[1] exp(-R2 TE)
  #           S_MT(TE,alpha) = th[2] exp(-R2 TE)
  #           S_PD(TE,alpha) = th[3] exp(-R2 TE)
  #      parameters th
  #   to      S_T1(TE,alpha) = p[1] exp(-R2 TE) sin(a1) (1-exp(-TR p[2]))/(1-cos(a1) exp(-TR p[2]))
  #           S_MT(TE,alpha) = p[1] exp(-R2 TE) p[3] sin(a2) (1-exp(-TR p[2]))/(1-cos(a1) exp(-TR p[2]))
  #           S_PD(TE,alpha) = p[1] exp(-R2 TE) sin(a2) (1-exp(-TR p[2]))/(1-cos(a2) exp(-TR p[2]))
  #      parameters  p
  #
  sa1th3 <- sin(a1)*th[3]
  sa2th1 <- sin(a2)*th[1]
  ca1 <- cos(a1)
  ca2 <- cos(a2)
  emTRR1 <- (sa2th3-sa1th1)/(ca1*sa2th3-ca2*sa1th1)
  R1 <-   -log(emTRR1)/TR
  S0 <- th[1]/sin(a1)*(1-ca1*emTRR1)/(1-emTRR1)
  M <- th[2]/th[3]
}

#   R1conf <- function(theta,si2,aT1,aPD,TR=1,df=NULL,alpha=0.05){
ESTATICS.confidence.old <- function(theta,si2,aT1,aPD,TR=1,df=NULL,alpha=0.05){
  ##
  ##   construct confidence region for parameter E1
  ##   parameters are expected to be named as "ST1", "SPD" and "R2star",
  ##          parameter "SMT" is not used
  ##  si2 - inverse covariance matrix of parameters
  ##
  fc1 <- sin(aT1)/sin(aPD)
  fc2 <- cos(aT1)
  fc3 <- fc1*cos(aPD)
  th <- theta[c("ST1","SPD")]
  E1 <- (th[1]-th[2]*fc1)/(th[1]*fc2-th[2]*fc3)
  R1 <- -log(E1)/TR
  R2star <- theta["R2star"]
  if(any(theta==0)){
    ## no interior solution, unable to provide confidence regions
    return(list(E1=E1,CIE1=c(NA,NA),R1=R1,CIR1=c(NA,NA),R2star=R2star,CIR2star=c(NA,NA)))
  }
  th <- theta[c("ST1","SPD")]
  Amat <- si2[c("ST1","SPD"),c("ST1","SPD")][c(1,2,4)]
  #    qnsq <- qnorm(1-alpha/2)^2
  qnsq <- if(is.null(df)) qchisq(1-alpha,2) else 2*qf(1-alpha,2,df)
  th2ofth1 <- function(th1,th,Amat,qnsq){
    th1diff <- (th1-th[1])
    p <- th[2]- th1diff*Amat[2]/Amat[3]
    q <- (qnsq - th1diff^2*Amat[1]+2*th1diff*th[2]*Amat[2])/Amat[3]-th[2]^2
    D <- p^2+q
    th2 <- if(D>=0) c(p-sqrt(D),p+sqrt(D)) else c(NA,NA)
    th2
  }

  ## now search for min/max of (1-fc1*th2(th1)/th1)/(fc2-fc3*th2(th1)/th1)
  e1 <- th[2]+th[1]*Amat[2]/Amat[3]
  e2 <- qnsq/Amat[3]
  e3 <- (Amat[1]*Amat[3]-Amat[2]^2)/Amat[3]^2
  e4 <- th[1]
  e0 <- e3^2*e4^2+e1^2*e3
  q <- (e3^2*e4^4+e2^2-2*e3*e4^2*e2-e1^2*e2+e1^2*e3*e4^2)/e0
  p <- -(e2*e3*e4-e3^2*e4^3-e1^2*e3*e4)/e0
  D <- p^2-q
  th1 <- if(D>=0) c(p-sqrt(D),p+sqrt(D)) else c(NA,NA)
  th2 <- th2ofth1(th1[2],th,Amat,qnsq)
  CIE1 <- (th1-th2*fc1)/(th1*fc2-th2*fc3)
  CIR1 <- -log(CIE1)/TR
  qqn <- if(is.null(df)) qnorm(1-alpha/2) else qt(1-alpha/2,df)
  CIR2star <- R2star + c(-qqn,qqn)/sqrt(si2["R2star","R2star"])
  list(E1=E1,CIE1=sort(CIE1),R1=R1,CIR1=sort(CIR1),R2star=R2star,CIR2star=CIR2star)
}

ESTATICS.confidence <- function(theta,si2,aT1,aPD,TR=1,df=NULL,alpha=0.05){
  ##
  ##   construct confidence region for parameter E1
  ##   parameters are expected to be named as "ST1", "SPD" and "R2star",
  ##          parameter "SMT" is not used
  ##  si2 - inverse covariance matrix of parameters
  ##
  fc1 <- sin(aT1)/sin(aPD)
  fc2 <- cos(aT1)
  fc3 <- fc1*cos(aPD)
  th <- theta[c("ST1","SPD")]
  E1 <- (th[1]-th[2]*fc1)/(th[1]*fc2-th[2]*fc3)
  R1 <- -log(E1)/TR
  R2star <- theta["R2star"]
  if(any(theta==0)){
    ## no interior solution, unable to provide confidence regions
    return(list(E1=E1,CIE1=c(NA,NA),R1=R1,CIR1=c(NA,NA),R2star=R2star,CIR2star=c(NA,NA)))
  }
  th <- theta[c("ST1","SPD")]
  Amat <- si2[c("ST1","SPD"),c("ST1","SPD")][c(1,2,4)]
  #    qnsq <- qnorm(1-alpha/2)^2
  qnsq <- if(is.null(df)) qchisq(1-alpha,2) else 2*qf(1-alpha,2,df)
  th2ofth1 <- function(th1,th,Amat,qnsq){
    th1diff <- (th1-th[1])
    p <- th[2]- th1diff*Amat[2]/Amat[3]
    q <- (qnsq - th1diff^2*Amat[1]+2*th1diff*th[2]*Amat[2])/Amat[3]-th[2]^2
    D <- p^2+q
    th2 <- if(D>=0) c(p-sqrt(D),p+sqrt(D)) else c(NA,NA)
    th2
  }

  ## now search for min/max of (1-fc1*th2(th1)/th1)/(fc2-fc3*th2(th1)/th1)
  e1 <- th[2]+th[1]*Amat[2]/Amat[3]
  e2 <- qnsq/Amat[3]
  e3 <- (Amat[1]*Amat[3]-Amat[2]^2)/Amat[3]^2
  e4 <- th[1]
  e0 <- e3^2*e4^2+e1^2*e3
  q <- (e3^2*e4^4+e2^2-2*e3*e4^2*e2-e1^2*e2+e1^2*e3*e4^2)/e0
  p <- -(e2*e3*e4-e3^2*e4^3-e1^2*e3*e4)/e0
  D <- p^2-q
  th1 <- if(D>=0) c(p-sqrt(D),p+sqrt(D)) else c(NA,NA)
  th2 <- th2ofth1(th1[2],th,Amat,qnsq)
  CIE1 <- (th1-th2*fc1)/(th1*fc2-th2*fc3)
  CIR1 <- -log(CIE1)/TR
  qqn <- if(is.null(df)) qnorm(1-alpha/2) else qt(1-alpha/2,df)
  CIR2star <- R2star + c(-qqn,qqn)/sqrt(si2["R2star","R2star"])
  list(E1=E1,CIE1=sort(CIE1),R1=R1,CIR1=sort(CIR1),R2star=R2star,CIR2star=CIR2star)
}

initth <- function(mpmdata){
  ## get initial estimates for ESTATICS parameters
  mask <- mpmdata$mask
  nvox <- prod(mpmdata$sdim)
  TE <- mpmdata$TE
  if(mpmdata$model==2){
     nT1 <- length(mpmdata$t1Files)
     indT1 <- c(1,nT1)
     nMT <- length(mpmdata$mtFiles)
     indMT <- nT1+c(1,nMT)
     nPD <- length(mpmdata$pdFiles)
     indPD <- nT1+nMT+c(1,nPD)
     th <- matrix(0,4,nvox)
     T1 <- matrix(mpmdata$ddata[indT1,,,],2,nvox)[,mask]
     MT <- matrix(mpmdata$ddata[indMT,,,],2,nvox)[,mask]
     PD <- matrix(mpmdata$ddata[indPD,,,],2,nvox)[,mask]
     R2star <- pmax(0,((log(T1[1,])-log(T1[2,]))/diff(TE[indT1])+
             (log(MT[1,])-log(MT[2,]))/diff(TE[indMT])+
             (log(PD[1,])-log(PD[2,]))/diff(TE[indPD]))/3)
     th[1,mask] <- T1[1,]*exp(R2star*TE[1])
     th[2,mask] <- MT[1,]*exp(R2star*TE[nT1+1])
     th[3,mask] <- PD[1,]*exp(R2star*TE[nT1+nMT+1])
     th[4,mask] <- R2star
     dim(th) <- c(4,mpmdata$sdim)
  }
  if(mpmdata$model==1){
     nT1 <- length(mpmdata$t1Files)
     indT1 <- c(1,nT1)
     n2 <- length(mpmdata$pdFiles)
     ind2 <- c(nT1+1,mpmdata$nFiles)
     th <- matrix(0,3,nvox)
     T1 <- matrix(mpmdata$ddata[indT1,,,],2,nvox)[,mask]
     S <- matrix(mpmdata$ddata[ind2,,,],2,nvox)[,mask]
     R2star <- pmax(0,((log(T1[1,])-log(T1[2,]))/diff(TE[indT1])+
             (log(S[1,])-log(S[2,]))/diff(TE[ind2]))/2)
     th[1,mask] <- T1[1,]*exp(R2star*TE[1])
     th[2,mask] <- S[1,]*exp(R2star*TE[nT1+1])
     th[3,mask] <- R2star
     dim(th) <- c(3,mpmdata$sdim)
  }
  if(mpmdata$model==0){
     nT1 <- mpmdata$nFiles
     indT1 <- c(1,nT1)
     th <- matrix(0,2,nvox)
     T1 <- matrix(mpmdata$ddata[indT1,,,],2,nvox)[,mask]
     R2star <- pmax(0,(log(T1[1,])-log(T1[2,]))/diff(TE[indT1]))
     th[1,mask] <- T1[1,]*exp(R2star*TE[1])
     th[2,mask] <- R2star
     dim(th) <- c(2,mpmdata$sdim)
  }
  th
}
