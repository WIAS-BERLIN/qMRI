\name{estimateIRsolidfixed}
\alias{estimateIRsolidfixed}
\title{
Estimate mixture parameter in Inversion Recovery MRI experiments mixture model for non-fluid voxel
}
\description{
Reestimate proportion of fluid with Sx and Rx fixed after smoothing.
}
\usage{
estimateIRsolidfixed(IRmixedobj, TEScale = 100, dataScale = 1000,
verbose = TRUE, lower = c(0), upper = c(0.95))
}
\arguments{
  \item{IRmixedobj}{ Object of class \code{"IRmixed"} as generated by function \code{\link{smoothIRSolid}} or \code{\link{estimateIRsolid}}.
}

  \item{TEScale}{
Internal scale factor for Echo Times. This influences parameter scales in numerical calculations.
}
  \item{dataScale}{
Internal scale factor for MR signals. This influences parameter scales in numerical calculations.
}
  \item{verbose}{
Logical. Provide some runtime diagnostics.
}
  \item{lower}{
   lower bound for fx (fluid proportion)
}
  \item{upper}{
   upper bound for fx (fluid proportion)
}
}
\value{
List of class \code{"IRmixed"} components
\item{IRdata }{4D array containing the IRMRI data, first dimension refers to inversion times}
\item{InvTimes }{vector of inversion times}
\item{segm}{segmentation codes, 1 for CSF, 2 for GM, 3 for WM, 0 for out of brain}
\item{sigma}{noise standard deviation, if not specified estimated fron CSF areas in image with largest inversion time}
\item{L}{effective number of coils}
\item{fx}{Array of fluid proportions}
\item{Sx}{Array of maximal signals}
\item{Rx}{Array of relaxation rates}
\item{Sf}{Global estimate of maximal fluid signal}
\item{Rf}{Global estimate of fluid relaxation rate}
\item{ICovx}{Covariance matrix of estimates \code{fx}, \code{Sx} and \code{Rx}.}
\item{sigma}{Array of provided or estimated noise standard deviations }
\item{Convx}{Array of convergence indicators}
\item{rsdx}{Residual standard deviations}
\item{method}{\code{"NLS"} for nonlinear regression or \code{"QL"} for quasi likelihood.}
\item{varest}{Method used for variance estimation}
The arrays contain entries for all voxel with \code{segments\%in\%1:3}.
}
\author{
Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{
 \code{\link{estimateIRfluid}}, \code{\link{estimateIRsolid}}, \code{\link{estimateIR}},\code{\link{smoothIRSolid}},
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (IRdata, InvTimes, segments, Sfluid, Rfluid, Ssolid, 
    Rsolid, TEScale = 100, dataScale = 1000, method = c("NLR", 
        "QL"), sigma = NULL, L = 1, maxR2star = 50, varest = c("RSS", 
        "data"), verbose = TRUE, lower = c(0.05), upper = c(0.95)) 
{
    mask <- segments > 1
    nvoxel <- sum(mask)
    ntimes <- length(InvTimes)
    InvTimes[InvTimes == Inf] <- 50 * max(InvTimes[InvTimes != 
        Inf])
    dimdata <- dim(IRdata)
    if (dimdata[1] != ntimes) 
        stop("estimateIRsolid: incompatible length of InvTimes")
    if (any(dimdata[-1] != dim(mask))) 
        stop("estimateIRsolid: incompatible dimension of segments")
    InvTimesScaled <- InvTimes/TEScale
    npar <- 1
    fx <- rsdx <- array(0, dim(mask))
    ICovx <- array(0, prod(dim(mask)))
    Convx <- array(0, dim(mask))
    fx[segments == 1] <- 1
    Rx <- Rsolid
    Sx <- Ssolid
    Convx[segments == 1] <- 0
    ICovx[segments == 1] <- 1e+20
    isConv <- array(FALSE, nvoxel)
    isThresh <- array(FALSE, nvoxel)
    modelCoeff <- numeric(nvoxel)
    invCov <- numeric(nvoxel)
    rsigma <- numeric(nvoxel)
    if (method == "QL") {
        if (is.null(sigma)) {
            method <- "NLR"
            warning("estimateIRsolid: method QL needs sigma estimated from fluid or supplied")
        }
        sig <- sigma/dataScale
        CL <- sig * sqrt(pi/2) * gamma(L + 0.5)/gamma(L)/gamma(1.5)
    }
    dim(IRdata) <- c(dimdata[1], prod(dim(segments)))
    IRdataSolid <- IRdata[, mask]
    Rsm <- Rsolid[mask]
    Ssm <- Ssolid[mask]
    thetas <- rep(0.1, nvoxel)
    if (verbose) {
        cat("Start estimation in", nvoxel, "voxel at", format(Sys.time()), 
            "\n")
        pb <- txtProgressBar(0, nvoxel, style = 3)
    }
    for (xyz in 1:nvoxel) {
        ivec <- IRdataSolid[, xyz]/dataScale
        th <- thetas[, xyz]
        Rs <- Rsm[xyz]
        Ss <- Ssm[xyz]
        res <- if (method == "NLR") 
            try(nls(ivec ~ IRmix2fix(par, ITS, Sf, Ss, Rf, Rs), 
                data = list(ITS = InvTimesScaled, Sf = Sfluid, 
                  Ss = Ss, Rf = Rfluid, Rs = Rs), start = list(par = th), 
                control = list(maxiter = 200, warnOnly = TRUE)), 
                silent = TRUE)
        else try(nls(ivec ~ IRmix2fixQL(par, ITS, Sf, Ss, Rf, 
            Rs, CL, sig, L), data = list(ITS = InvTimesScaled, 
            Sf = Sfluid, Ss = Ss, Rf = Rfluid, Rs = Rs, CL = CL, 
            sig = sig, L = L), start = list(par = th), control = list(maxiter = 200, 
            warnOnly = TRUE)), silent = TRUE)
        if (class(res) == "try-error") {
            th <- pmin(upper, pmax(lower, th))
            res <- if (method == "NLR") 
                try(nls(ivec ~ IRmix2fix(par, ITS, Sf, Ss, Rf, 
                  Rs), data = list(ITS = InvTimesScaled, Sf = Sfluid, 
                  Ss = Ss, Rf = Rfluid, Rs = Rs), start = list(par = th), 
                  algorithm = "port", control = list(maxiter = 200, 
                    warnOnly = TRUE), lower = lower, upper = upper), 
                  silent = TRUE)
            else try(nls(ivec ~ IRmix2fixQL(par, ITS, Sf, Ss, 
                Rf, Rs, CL, sig, L), data = list(ITS = InvTimesScaled, 
                Sf = Sfluid, Ss = Ss, Rf = Rfluid, Rs = Rs, CL = CL, 
                sig = sig, L = L), start = list(par = th), algorithm = "port", 
                control = list(maxiter = 200, warnOnly = TRUE), 
                lower = lower, upper = upper), silent = TRUE)
        }
        if (class(res) != "try-error") {
            sres <- if (varest == "RSS") 
                getnlspars(res)
            else getnlspars2(res, shat[, xyz], sind)
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
    list(fx = fx, Rx = Rx, Sx = Sx, Sf = Sfluid, Rf = Rfluid, 
        ICovx = ICovx, Convx = Convx, sigma = sigma, rsdx = rsdx)
  }
}
 \keyword{ models}
 \keyword{ regression }
\concept{ IRMRI }
