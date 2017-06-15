\name{estimateESTATICS}
\alias{estimateESTATICS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate parameters in the ESTATICS model using nonlinear least squares regression.
}
\description{
Evaluation of the 
ESTATICS model (Weisskopf (2013) using nonlinear least squares regression 
and the second a quasi-likelihood approach assuming a noncentral chi- or a Rician distribuion
for the data. The latter function should be preferred in case of low SNR (high resolution) 
data to avoid biased parameter estimates. Function \code{\link{estimateESTATICSQL}} requires
a specification of the scale parameter sigma of the data distribution. For the case of 
the Rician distribution paramter maps of sigma can be estimated from magnitude and phase images 
(complex polar coordinates) using function \code{\link{estimateSigma}}.
}
\usage{
estimateESTATICS(mpmdata, TEScale = 100, dataScale = 1000, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mpmdata}{
  Object containing the measured MPM data.
}
  \item{TEScale}{
  scale factor for TE (used for improved numerical stability)
}
  \item{dataScale}{
  scale factor for image intensities (used for improved numerical stability)
}
  \item{verbose}{
  logical: Monitor process.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 list with components
                 \item{modelCoeff }{Estimated parameter maps}
                 \item{invCov }{ map of inverse covariance matrices}
                 \item{isConv }{ convergence indicator map}
                 \item{sdim }{ image dimension}
                 \item{nFiles }{ number of images}
                 \item{t1Files }{ vector of T1 filenames}
                 \item{pdFiles }{ vector of PD filenames}
                 \item{mtFiles }{ vector of MT filenames}
                 \item{model }{ model used (depends on specification of MT files)}
                 \item{maskFile }{ filename of brain mask}
                 \item{mask }{ brain mask}
                 \item{TR }{ TR values}
                 \item{TE }{ TE values}
                 \item{FA }{ Flip angles (FA)}
                 \item{TEScale }{ TEScale}
                 \item{dataScale }{dataScale}
}
\references{
Weiskopf, N.; Suckling, J.; Williams, G.; Correia, M. M.; Inkster, B.; Tait, R.; Ooi, C.; Bullmore, E. T. & Lutti, A. Quantitative multi-parameter mapping of R1, PD(*), MT, and R2(*) at 3T: a multi-center validation. Front Neurosci, Wellcome Trust Centre for Neuroimaging, UCL Institute of Neurology, University College London London, UK., 2013, 7, 95
}
\author{
Karsten Tabelow \email{karsten.tabelow@wias-berlin.de}, Joerg Polzehl \email{joerg.polzehl@wias-berlin.de}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{estimateQI}},\code{\link{estimateQIconf}},\code{\link{estimateESTATICSQL}},\code{\link{smoothESTATICS}}
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
