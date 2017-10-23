\name{estimateQIconf}
\alias{estimateQIconf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get confidence intervall for QI parameters
}
\description{
Quantitaive imaging parameters are calculated from the estimated parameters 
in the ESTATICS model. This involves a correction for magnetic field inhomogeneities
if the information is provided in argument \code{b1File} and use of a second 
of a second recovery delay TR2 in case of Dual-Exitation FLASH measurements
(Helms 2008).
}
\usage{
estimateQIconf(mpmESTATICSmodel, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mpmESTATICSmodel}{
%%     ~~Describe \code{mpmESTATICSmodel} here~~
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
function (mpmESTATICSmodel, verbose = TRUE) 
{
    if (mpmESTATICSmodel$model != 2) 
        stop("only full model implemented!")
    dimnames(mpmESTATICSmodel$modelCoeff) <- list(c("ST1", "SMT", 
        "SPD", "R2star"), NULL, NULL, NULL)
    dimnames(mpmESTATICSmodel$invCov) <- list(c("ST1", "SMT", 
        "SPD", "R2star"), c("ST1", "SMT", "SPD", "R2star"), NULL, 
        NULL, NULL)
    R1 <- array(0, mpmESTATICSmodel$sdim)
    CIR1 <- array(0, c(2, mpmESTATICSmodel$sdim))
    R2 <- array(0, mpmESTATICSmodel$sdim)
    CIR2 <- array(0, c(2, mpmESTATICSmodel$sdim))
    for (z in 1:mpmESTATICSmodel$sdim[3]) {
        for (y in 1:mpmESTATICSmodel$sdim[2]) {
            for (x in 1:mpmESTATICSmodel$sdim[1]) {
                if (mpmESTATICSmodel$mask[x, y, z]) {
                  if (is.null(mpmESTATICSmodel$bi)) {
                    zz <- ESTATICS.confidence(mpmESTATICSmodel$modelCoeff[, 
                      x, y, z], mpmESTATICSmodel$invCov[, , x, 
                      y, z], mpmESTATICSmodel$FA[1] * pi/180, 
                      mpmESTATICSmodel$FA[length(mpmESTATICSmodel$t1Files) + 
                        length(mpmESTATICSmodel$mtFiles) + 1] * 
                        pi/180, mpmESTATICSmodel$TR[1], df = mpmESTATICSmodel$nFiles - 
                        4, 0.05)
                  }
                  else {
                    zz <- ESTATICS.confidence(mpmESTATICSmodel$modelCoeff[, 
                      x, y, z], mpmESTATICSmodel$invCov[, , x, 
                      y, z] * mpmESTATICSmodel$bi[x, y, z], mpmESTATICSmodel$FA[1] * 
                      pi/180, mpmESTATICSmodel$FA[length(mpmESTATICSmodel$t1Files) + 
                      length(mpmESTATICSmodel$mtFiles) + 1] * 
                      pi/180, mpmESTATICSmodel$TR[1], df = mpmESTATICSmodel$nFiles - 
                      4, 0.05)
                  }
                  R2[x, y, z] <- zz$R2star
                  R1[x, y, z] <- zz$R1
                  CIR1[, x, y, z] <- zz$CIR1
                  CIR2[, x, y, z] <- zz$CIR2star
                }
            }
        }
        if (verbose) 
            cat(z, format(Sys.time()), "\n")
    }
    invisible(list(R1 = R1, CIR1 = CIR1, R2star = R2, CIR2star = CIR2, 
        model = mpmESTATICSmodel$model, t1Files = mpmESTATICSmodel$t1Files, 
        mtFiles = mpmESTATICSmodel$mtFiles, pdFiles = mpmESTATICSmodel$pdFiles, 
        mask = mpmESTATICSmodel$mask))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
