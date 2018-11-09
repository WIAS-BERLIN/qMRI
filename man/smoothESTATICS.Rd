\name{smoothESTATICS}
\alias{smoothESTATICS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Structural adaptive moothing of parameter maps in the ESTATICS model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
smoothESTATICS(mpmESTATICSModel, mpmData = NULL, kstar = 16, alpha = 0.0025,
  patchsize=0, wghts = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mpmESTATICSModel}{
%%     ~~Describe \code{mpmESTATICSModel} here~~
}
  \item{mpmData}{
%%     ~~Describe \code{mpmData} here~~
}
  \item{kstar}{
%%     ~~Describe \code{kstar} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{patchsize}{
%%     ~~Describe \code{patchsize} here~~
}
  \item{wghts}{
%%     ~~Describe \code{wghts} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (mpmESTATICSModel, mpmData = NULL, kstar = 16, alpha = 0.05,
    wghts = NULL, verbose = TRUE)
{
    nv <- if (mpmESTATICSModel$model == 2)
        4
    else 3
    lambda <- nv * qf(1 - alpha, nv, mpmESTATICSModel$nFiles -
        nv)
    if (is.null(wghts))
        wghts <- c(1, 1, 1)
    wghts <- wghts[1]/wghts[2:3]
    n1 <- mpmESTATICSModel$sdim[1]
    n2 <- mpmESTATICSModel$sdim[2]
    n3 <- mpmESTATICSModel$sdim[3]
    n <- n1 * n2 * n3
    zobj <- list(bi = rep(1, n), theta = mpmESTATICSModel$modelCoeff)
    bi <- zobj$bi
    mc.cores <- setCores(, reprt = FALSE)
    hmax <- 1.25^(kstar/3)
    if (verbose) {
        mae <- NULL
        protocol <- matrix("", kstar, 1, dimnames = list(paste("step",
            1:kstar), "protocol"))
    }
    k <- 1
    if (verbose)
        pb <- txtProgressBar(min = 0, max = kstar, style = 3)
    if (!is.null(mpmData)) {
        if (length(dim(mpmData$ddata)) != 4 | any(mpmESTATICSModel$sdim !=
            dim(mpmData$ddata)[-1]))
            stop("incompatible dimensions of model parameters and original data")
        kstar <- kstar - 1
        smoothData <- TRUE
    }
    else {
        smoothData <- FALSE
    }
    while (k <= kstar) {
        hakt <- gethani(1, 1.25 * hmax, 2, 1.25^k, wghts, 1e-04)
        dlw <- (2 * trunc(hakt/c(1, wghts)) + 1)[1:3]
        zobj <- .Fortran("vaws", as.double(mpmESTATICSModel$modelCoeff),
            as.logical(mpmESTATICSModel$mask), as.integer(nv),
            as.integer(n1), as.integer(n2), as.integer(n3), hakt = as.double(hakt),
            as.double(lambda), as.double(zobj$theta), as.double(mpmESTATICSModel$invCov),
            bi = as.double(zobj$bi), theta = double(nv * n),
            as.integer(mc.cores), double(prod(dlw)), as.double(wghts),
            double(nv * mc.cores), PACKAGE = "qMRI")[c("bi",
            "theta", "hakt")]
        bi <- zobj$bi <- pmax(bi, zobj$bi)
        if (verbose) {
            protocol[k] <- paste("bandwidth: ", signif(hakt,
                3), "MSE: ", signif(mean((zobj$theta - mpmESTATICSModel$modelCoeff)^2),
                3), "MAE: ", m1 <- signif(mean(abs(zobj$theta -
                mpmESTATICSModel$modelCoeff)), 3), "mean(bi):",
                signif(mean(zobj$bi), 3))
            mae <- c(mae, m1)
            setTxtProgressBar(pb, k)
        }
        k <- k + 1
        gc()
    }
    if (smoothData) {
        nve <- mpmESTATICSModel$nFiles
        hakt <- gethani(1, 1.25 * hmax, 2, 1.25^k, wghts, 1e-04)
        dlw <- (2 * trunc(hakt/c(1, wghts)) + 1)[1:3]
        zobj <- .Fortran("vawsext", as.double(mpmESTATICSModel$modelCoeff),
            as.logical(mpmESTATICSModel$mask), as.integer(nv),
            as.integer(n1), as.integer(n2), as.integer(n3), as.double(mpmData$ddata),
            as.integer(nve), hakt = as.double(hakt), as.double(lambda),
            as.double(zobj$theta), as.double(mpmESTATICSModel$invCov),
            bi = as.double(zobj$bi), theta = double(nv * n),
            thext = double(nve * n), as.integer(mc.cores), double(prod(dlw)),
            as.double(wghts), double(nv * mc.cores), double(nve *
                mc.cores), PACKAGE = "qMRI")[c("bi", "theta",
            "thext", "hakt")]
        dim(zobj$thext) <- c(mpmESTATICSModel$nFiles, mpmESTATICSModel$sdim)
        mpmDataSmoothed <- zobj$thext
        zobj$bi <- pmax(bi, zobj$bi)
        if (verbose) {
            protocol[k] <- paste("bandwidth: ", signif(hakt,
                3), "MSE: ", signif(mean((zobj$theta - mpmESTATICSModel$modelCoeff)^2),
                3), "MAE: ", m1 <- signif(mean(abs(zobj$theta -
                mpmESTATICSModel$modelCoeff)), 3), "mean(bi):",
                signif(mean(zobj$bi), 3))
            mae <- c(mae, m1)
            setTxtProgressBar(pb, kstar)
        }
    }
    else {
        mpmDataSmoothed <- NULL
    }
    if (verbose)
        close(pb)
    if (verbose)
        print(protocol)
    dim(zobj$theta) <- c(nv, n1, n2, n3)
    dim(zobj$bi) <- c(n1, n2, n3)
    invisible(list(modelCoeff = zobj$theta, invCov = mpmESTATICSModel$invCov,
        isConv = mpmESTATICSModel$isConv, bi = zobj$bi, smoothPar = c(lambda,
            hakt, alpha), smoothedData = mpmDataSmoothed, sdim = mpmESTATICSModel$sdim,
        nFiles = mpmESTATICSModel$nFiles, t1Files = mpmESTATICSModel$t1Files,
        pdFiles = mpmESTATICSModel$pdFiles, mtFiles = mpmESTATICSModel$mtFiles,
        model = mpmESTATICSModel$model, maskFile = mpmESTATICSModel$maskFile,
        mask = mpmESTATICSModel$mask, TR = mpmESTATICSModel$TR,
        TE = mpmESTATICSModel$TE, FA = mpmESTATICSModel$FA, TEScale = mpmESTATICSModel$TEScale,
        dataScale = mpmESTATICSModel$dataScale))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
