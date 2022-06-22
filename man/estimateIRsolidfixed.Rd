\name{estimateIRsolidfixed}
\alias{estimateIRsolidfixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate mixture parameter in Inversion Recovery MRI experiments mixture model for non-fluid voxel
}
\description{
Reestimate proportion of fluid with Sx and Rx fixed after smoothing.
}
\usage{
estimateIRsolidfixed(IRdata, InvTimes, segments, Sfluid, Rfluid,
Ssolid, Rsolid, TEScale = 100, dataScale = 1000, method = c("NLR", "QL"),
sigma = NULL, L = 1, varest = c("RSS", "data"),
verbose = TRUE, lower = c(0), upper = c(0.95))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IRdata}{ 4D array of IRMRI signals. First dimension corresponds to Inversion times (InvTime).
}
  \item{InvTimes}{
  Inversion Times used in the experiments
}d
  \item{segments}{
  Results of tissue segmentation using, e.g., function \code{fast} from package \code{fslr}. The array is assumed to contain values 1 for CSF voxel, 2 for GM and 3 for WM.
}
  \item{Sfluid}{
  maximal signal for fluid compartment
}
  \item{Rfluid}{
  Recovery rate for fluid compartments
}
  \item{Ssolid}{
  Array of signals for pure solid/fluid
}
  \item{Rsolid}{
  Array of Recovery rates
}
  \item{TEScale}{
Internal scale factor for Echo Times. This influences parameter scales in numerical calculations.
}
  \item{dataScale}{
Internal scale factor for MR signals. This influences parameter scales in numerical calculations.
}
  \item{method}{
  Either \code{"NLS"} for nonlinear least squares (ignores Rician bias) or \code{"QL"} for Quasi-Likelihood. The second option is more accurate but requires additional information and is computationally more expensive.
}
  \item{sigma}{
Noise standard deviation of (unprocessed) MRI signals. This needs to be accessed from minimally processed data or directly calculated from the scanner protocol. Only needed in case of \code{method=="QL"}.
}d
  \item{L}{
Number of effective coils. Depends on the reconstruction method on multi-coil MR systems. Should be 1 in case of SENSE ceconstructions. In case of GRAPPA this varies spatially and needs to be calculated.
}
  \item{varest}{
  Method to, in case of \code{method="QR"}, estimate \code{sigma}if not provided. Either from residual sums of squares (\code{"RSS"}) or MR signals (\code{"data"}) using function \code{varest} gtom package {\code{aws}}. Only to be used in case that no image registration was needed as preprocessing.
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
List of components
\item{fx}{Array of fluid proportions}
\item{Sx}{Array of maximal signals}
\item{Rx}{Array of relaxation rates}
\item{Sf}{Global estimate of maximal fluid signal}
\item{Rf}{Global estimate of fluid relaxation rate}
\item{ICovx}{Covariance matrix of estimates \code{fx}, \code{Sx} and \code{Rx}.}
\item{sigma}{Array of provided or estimated noise standard deviations }
\item{Convx}{Array of convergence indicators}
\item{rsdx}{Residual standard deviations}
The arrays contain entries for all voxel with \code{segments\%in\%1:3}.
}
\author{
Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
J\"org Polzehl \email{polzehl@wias-berlin.de}
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{estimateIRfluid}}, \code{\link{estimateIRsolid}}, \code{\link{estimateIR}},\code{\link{smoothIRSolid}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
