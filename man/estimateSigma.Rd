\name{estimateSigma}
\alias{estimateSigma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate scale parameter maps for Rician distributed data from 3D magnitude and phase images.
}
\description{
Magnitude and phase images contain the measured complex image information after image reconstruction using e.g. SENSE
in polar coordinates. This information is transformed in images containing the corresponding real and imaginary images.
A vectorized adaptive weights algorithm (Polzehl and Spokoiny, 2006) is used to estimate the complex image intensities.
In each step magnitude maps of the complex normalized residuals are obtained (these follow a Rician distribution) and smoothed using a spatial median filter with bandwidth \code{hsig} to obtain estimated maps of scale parameter sigma. The final magnitude maps of the complex normalized residuals, magnitude (intensity) maps and maps of sigma are returned.




}
\usage{
estimateSigma(magnitude, phase, mask, kstar = 20, kmin = 10, hsig = 5, lambda = 12, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{magnitude}{
  Either a file name corresponding to a 3D magnitude image in Nifti (.nii) format or a 3D array containing the
  image data.
}
  \item{phase}{
  Either a file name corresponding to a 3D phase image in Nifti (.nii) format or a 3D array containing the
  image data.
}
  \item{mask}{
  Either a logical 3D array describing a brain mask or \code{NULL} 
}
  \item{kstar}{
  Maximum number of steps in the adaptive smoothing algorithm used to estimate the spatial mean intensity map.
}
  \item{kmin}{
  Number of steps without adaptation (the first \code{kmin-1} steps are not carried out.
}
  \item{hsig}{
  bandwidth for spatial median filtering of sum of squared residuals
}
  \item{lambda}{
  adaptation parameter (chosen by a propagation condition), smaller values may lead to severe underestimation
  of the scale parameter due to adaptation to the data, significantly larger values may lead to
  positively biased estimates. Please do not change ...
}
  \item{verbose}{
  Logical: Monitor the progress.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class 'sigmaEstSENSE' with components
       \item{sigma}{3D scale parameter map}
       \item{theta}{Magnitude image estimated by adaptive smoothing}
       \item{sigmal}{Array of residual magnitudes}
       \item{mask}{Brain mask used} 
       \item{protocol}{Protocoll of iteration process}
       \item{args}{Arguments of call}
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Joerg Polzehl \email{joerg.polzehl@wias-berlin.de}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{estimateESTATICSQL}}, \code{\link{medianFilterSigma}} %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (magnitude, phase, mask, kstar = 20, kmin = 10, hsig = 5, 
    lambda = 12, verbose = TRUE) 
{
    args <- sys.call(-1)
    sdim <- dim(mask)
    if (!is.numeric(magnitude)) {
        if (verbose) 
            cat("reading Magnitude file ... ")
        R <- readNIfTI(magnitude, reorient = FALSE)
    }
    else {
        R <- magnitude
    }
    if (!is.numeric(phase)) {
        if (verbose) 
            cat("reading Phase file ... ")
        Ph <- readNIfTI(phase, reorient = FALSE)
    }
    else {
        Ph <- phase
    }
    ComplImg <- array(0, c(2, sdim))
    ComplImg[1, , , ] <- R * cos(Ph)
    ComplImg[2, , , ] <- R * sin(Ph)
    mc.cores <- setCores(, reprt = FALSE)
    n <- prod(sdim)
    lambda0 <- 1e+40
    sigma2 <- array(1e+10, sdim)
    k <- kmin
    hmax <- 1.25^(kstar/3)
    nwmd <- (2 * as.integer(hsig) + 1)^3
    parammd <- .Fortran("paramw3", as.double(hsig), as.double(c(1, 
        1)), ind = integer(3 * nwmd), w = double(nwmd), n = as.integer(nwmd), 
        PACKAGE = "qMRI")[c("ind", "w", "n")]
    nwmd <- parammd$n
    parammd$ind <- parammd$ind[1:(3 * nwmd)]
    dim(parammd$ind) <- c(3, nwmd)
    if (verbose) 
        pb <- txtProgressBar(min = 0, max = kstar - kmin + 1, 
            style = 3)
    bi <- array(1, sdim)
    zobj <- list(theta = ComplImg, bi = bi)
    if (verbose) {
        mae <- NULL
        protocol <- matrix("", kstar - kmin + 1, 1, dimnames = list(paste("step", 
            kmin:kstar), "protocol"))
    }
    while (k <= kstar) {
        hakt <- gethani(1, 1.25 * hmax, 2, 1.25^k, c(1, 1), 1e-04)
        dlw <- (2 * trunc(hakt/c(1, 1, 1)) + 1)[1:3]
        zobj <- .Fortran("vaws2", as.double(ComplImg), as.logical(mask), 
            as.integer(2), as.integer(sdim[1]), as.integer(sdim[2]), 
            as.integer(sdim[3]), hakt = as.double(hakt), as.double(lambda0), 
            as.double(zobj$theta), as.double(sigma2), bi = as.double(zobj$bi), 
            theta = double(2 * n), sigma2 = double(n), as.integer(mc.cores), 
            double(prod(dlw)), as.double(c(1, 1)), double(2 * 
                mc.cores), PACKAGE = "qMRI")[c("bi", "theta", 
            "hakt", "sigma2")]
        dim(zobj$sigma2) <- sdim
        sigma2 <- .Fortran("mediansm", as.double(zobj$sigma2), 
            as.logical(mask), as.integer(sdim[1]), as.integer(sdim[2]), 
            as.integer(sdim[3]), as.integer(parammd$ind), as.integer(nwmd), 
            double(nwmd * mc.cores), as.integer(mc.cores), sigma2n = double(n), 
            PACKAGE = "qMRI")$sigma2n/0.6931
        bi <- zobj$bi <- pmax(bi, zobj$bi)
        if (verbose) {
            protocol[k - kmin + 1, 1] <- paste("bandwidth: ", 
                signif(hakt, 3), "sigma: mean: ", signif(sqrt(mean(sigma2[mask])), 
                  3), "median: ", signif(sqrt(median(sigma2[mask])), 
                  3), "sd: ", signif(sd(sqrt(sigma2[mask])), 
                  3), "median(bi):", signif(median(zobj$bi[mask]), 
                  3), "max(bi):", signif(max(zobj$bi[mask]), 
                  3))
            setTxtProgressBar(pb, k - kmin + 1)
        }
        k <- k + 1
        lambda0 <- lambda
        gc()
    }
    dim(zobj$theta) <- c(2, sdim)
    z <- list(sigma = array(sqrt(sigma2), sdim), theta = array(sqrt(zobj$theta[1, 
        , , ]^2 + zobj$theta[2, , , ]^2), sdim), sigmal = array(sqrt(zobj$sigma2), 
        sdim), mask = mask, protocol = protocol, args = args)
    class(z) <- "sigmaEstSENSE"
    z
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
