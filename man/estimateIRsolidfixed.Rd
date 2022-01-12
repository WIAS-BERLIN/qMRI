\name{estimateIRsolidfixed}
\alias{estimateIRsolidfixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate mixture parameter in Inversion Recovery MRI experiments mixture model for non-fluid voxel
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimateIRsolidfixed(IRdata, InvTimes, segments, Sfluid, Rfluid, Ssolid, Rsolid, TEScale = 100, dataScale = 1000, method = c("NLR", "QL"), sigma = NULL, L = 1, maxR2star = 50, varest = c("RSS", "data"), verbose = TRUE, lower = c(0.05), upper = c(0.95))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IRdata}{
%%     ~~Describe \code{IRdata} here~~
}
  \item{InvTimes}{
%%     ~~Describe \code{InvTimes} here~~
}
  \item{segments}{
%%     ~~Describe \code{segments} here~~
}
  \item{Sfluid}{
%%     ~~Describe \code{Sfluid} here~~
}
  \item{Rfluid}{
%%     ~~Describe \code{Rfluid} here~~
}
  \item{Ssolid}{
%%     ~~Describe \code{Ssolid} here~~
}
  \item{Rsolid}{
%%     ~~Describe \code{Rsolid} here~~
}
  \item{TEScale}{
%%     ~~Describe \code{TEScale} here~~
}
  \item{dataScale}{
%%     ~~Describe \code{dataScale} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{maxR2star}{
%%     ~~Describe \code{maxR2star} here~~
}
  \item{varest}{
%%     ~~Describe \code{varest} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{lower}{
%%     ~~Describe \code{lower} here~~
}
  \item{upper}{
%%     ~~Describe \code{upper} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
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
