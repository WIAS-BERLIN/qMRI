\name{medianFilterSigma}
\alias{medianFilterSigma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use a spatial median filter on residual sum of squares maps in objects
produced by function \code{estimateSigma}.
}
\description{
Apply a spatial median filter with bandwidth \code{hsig} to an 3D array or to residual magnitude map in component \code{sigmal}
of an object of class 'sigmaEstSENSE' returned by function \code{\link{estimateSigma}}
}
\usage{
medianFilterSigma(obj, hsig = 10, mask = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
     3D array or object of class 'sigmaEstSENSE'
}
  \item{hsig}{
     bandwidth for median filter
}
  \item{mask}{
     brain mask
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An array or an object of class 'sigmaEstSENSE' with components
       \item{sigma}{3D scale parameter map (replaced by result)}
       \item{theta}{Magnitude image estimated by adaptive smoothing}
       \item{sigmal}{Array of residual magnitudes}
       \item{mask}{Brain mask used} 
       \item{protocol}{Protocoll of iteration process}
       \item{args}{Arguments of call}
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
\code{\link{estimateSigma}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (obj, hsig = 10, mask = NULL) 
{
    if (class(obj) == "sigmaEstSENSE") {
        sigma2 <- obj$sigmal^2
        mask <- obj$mask
    }
    else {
        sigma2 <- obj^2
    }
    sdim <- dim(sigma2)
    n <- prod(sdim)
    if (length(sdim) != 3) 
        stop("obj needs to be of class 'array' (3D) or 'sigmaEstSENSE'")
    if (is.null(mask)) 
        mask <- array(TRUE, sdim)
    if (any(dim(mask) != sdim)) 
        stop("dimensions do not coinside")
    nwmd <- (2 * as.integer(hsig) + 1)^3
    parammd <- .Fortran("paramw3", as.double(hsig), as.double(c(1, 
        1)), ind = integer(3 * nwmd), w = double(nwmd), n = as.integer(nwmd), 
        PACKAGE = "qMRI")[c("ind", "w", "n")]
    nwmd <- parammd$n
    parammd$ind <- parammd$ind[1:(3 * nwmd)]
    dim(parammd$ind) <- c(3, nwmd)
    mc.cores <- setCores(, reprt = FALSE)
    sigma2 <- .Fortran("mediansm", as.double(sigma2), as.logical(mask), 
        as.integer(sdim[1]), as.integer(sdim[2]), as.integer(sdim[3]), 
        as.integer(parammd$ind), as.integer(nwmd), double(nwmd * 
            mc.cores), as.integer(mc.cores), sigma2n = double(n), 
        PACKAGE = "qMRI")$sigma2n/0.6931
    dim(sigma2) <- sdim
    if (class(obj) == "sigmaEstSENSE") {
        obj$sigma2 <- sqrt(sigma2)
        obj$hsig <- hsig
    }
    else {
        obj <- sqrt(sigma2)
    }
    obj
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utility }
\keyword{ smoothing }% __ONLY ONE__ keyword per line
