\name{estimateESTATICS}
\alias{estimateESTATICS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate parameters in the ESTATICS model using nonlinear least squares regression.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimateESTATICS(mpmdata, TEScale = 100, dataScale = 1000, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mpmdata}{
%%     ~~Describe \code{mpmdata} here~~
}
  \item{TEScale}{
%%     ~~Describe \code{TEScale} here~~
}
  \item{dataScale}{
%%     ~~Describe \code{dataScale} here~~
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
function (mpmdata, TEScale = 100, dataScale = 1000, verbose = TRUE) 
{
    if (mpmdata$model == 2) {
        xmat <- matrix(0, mpmdata$nFiles, 4)
        xmat[1:length(mpmdata$t1Files), 1] <- 1
        xmat[(length(mpmdata$t1Files) + 1):(length(mpmdata$t1Files) + 
            length(mpmdata$mtFiles)), 2] <- 1
        xmat[(length(mpmdata$t1Files) + length(mpmdata$mtFiles) + 
            1):mpmdata$nFiles, 3] <- 1
        xmat[, 4] <- mpmdata$TE/TEScale
        zerovoxel <- as.logical(((apply(mpmdata$ddata[xmat[, 
            1] == 1, , , ], 2:4, sum) == 0) | (apply(mpmdata$ddata[xmat[, 
            2] == 1, , , ], 2:4, sum) == 0) | (apply(mpmdata$ddata[xmat[, 
            3] == 1, , , ], 2:4, sum) == 0)) * mpmsense$mask)
    }
    else {
        xmat <- matrix(0, mpmdata$nFiles, 3)
        xmat[1:length(mpmdata$t1Files), 1] <- 1
        xmat[(length(mpmdata$t1Files) + 1):mpmdata$nFiles, 2] <- 1
        xmat[, 3] <- mpmdata$TE/TEScale
        zerovoxel <- as.logical(((apply(mpmdata$ddata[xmat[, 
            1] == 1, , , ], 2:4, sum) == 0) | (apply(mpmdata$ddata[xmat[, 
            2] == 1, , , ], 2:4, sum) == 0)) * mpmsense$mask)
    }
    if (verbose) {
        cat("Design of the model:\n")
        print(xmat)
    }
    mpmdata$mask[zerovoxel] <- FALSE
    R2star <- 0.05 * TEScale
    indT1 <- order(mpmdata$TE[as.logical(xmat[, 1])])[1]
    if (mpmdata$model == 2) {
        indMT <- order(mpmdata$TE[as.logical(xmat[, 2])])[1] + 
            sum(xmat[, 1])
        indPD <- order(mpmdata$TE[as.logical(xmat[, 3])])[1] + 
            sum(xmat[, 1]) + sum(xmat[, 2])
        npar <- 4
    }
    else {
        indPD <- order(mpmdata$TE[as.logical(xmat[, 2])])[1] + 
            sum(xmat[, 1])
        npar <- 3
    }
    isConv <- array(FALSE, mpmdata$sdim)
    modelCoeff <- array(0, c(npar, mpmdata$sdim))
    invCov <- array(0, c(npar, npar, mpmdata$sdim))
    if (verbose) 
        Sys.time()
    for (z in 1:mpmdata$sdim[3]) {
        for (y in 1:mpmdata$sdim[2]) {
            for (x in 1:mpmdata$sdim[1]) {
                if (mpmdata$mask[x, y, z]) {
                  ivec <- mpmdata$ddata[, x, y, z]/dataScale
                  if (mpmdata$model == 2) {
                    th <- c(ivec[indT1] * exp(-xmat[indT1, 4] * 
                      R2star), ivec[indMT] * exp(-xmat[indMT, 
                      4] * R2star), ivec[indPD] * exp(-xmat[indPD, 
                      4] * R2star), R2star)
                    res <- try(nls(ivec ~ qflashpl(par, xmat), 
                      start = list(par = th), control = list(maxiter = 200, 
                        warnOnly = TRUE)))
                    if (class(res) == "try-error" || !res$convInfo$isConv || 
                      any(coefficients(res) < 0)) 
                      res <- nls(ivec ~ qflashpl(par, xmat), 
                        start = list(par = th), algorithm = "port", 
                        control = list(warnOnly = TRUE, printEval = TRUE), 
                        lower = rep(0, 4))
                  }
                  else {
                    th <- c(ivec[indT1] * exp(-xmat[indT1, 3] * 
                      R2star), ivec[indPD] * exp(-xmat[indPD, 
                      3] * R2star), R2star)
                    res <- try(nls(ivec ~ qflashpl2(par, xmat), 
                      start = list(par = th), control = list(maxiter = 200, 
                        warnOnly = TRUE)))
                    if (class(res) == "try-error" || !res$convInfo$isConv || 
                      any(coefficients(res) < 0)) 
                      res <- nls(ivec ~ qflashpl2(par, xmat), 
                        start = list(par = th), algorithm = "port", 
                        control = list(warnOnly = TRUE, printEval = TRUE), 
                        lower = rep(0, 3))
                  }
                  sres <- summary(res)
                  isConv[x, y, z] <- res$convInfo$isConv
                  modelCoeff[, x, y, z] <- sres$coefficients[, 
                    1]
                  if (sres$sigma != 0) {
                    invCov[, , x, y, z] <- solve(sres$cov.unscaled)/sres$sigma^2
                  }
                }
            }
        }
        if (verbose) 
            cat(z, format(Sys.time()), "\n")
    }
    if (verbose) 
        Sys.time()
    invisible(list(modelCoeff = modelCoeff, invCov = invCov, 
        isConv = isConv, sdim = mpmdata$sdim, nFiles = mpmdata$nFiles, 
        t1Files = mpmdata$t1Files, pdFiles = mpmdata$pdFiles, 
        mtFiles = mpmdata$mtFiles, model = mpmdata$model, maskFile = mpmdata$maskFile, 
        mask = mpmdata$mask, TR = mpmdata$TR, TE = mpmdata$TE, 
        FA = mpmdata$FA, TEScale = TEScale, dataScale = dataScale))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
