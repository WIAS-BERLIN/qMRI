\name{calculateQI}
\alias{calculateQI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate quantitative imaging parameters from estimates in the ESTATICS model.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculateQI(mpmESTATICSModel, b1File = NULL, TR2 = 0, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mpmESTATICSModel}{
%%     ~~Describe \code{mpmESTATICSModel} here~~
}
  \item{b1File}{
%%     ~~Describe \code{b1File} here~~
}
  \item{TR2}{
%%     ~~Describe \code{TR2} here~~
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
function (mpmESTATICSModel, b1File = NULL, TR2 = 0, verbose = TRUE) 
{
    if (!is.null(b1File)) {
        if (verbose) 
            cat("reading B1 correction file from", b1File, "\n")
        b1Map <- readNIfTI(b1File, reorient = FALSE)/100
        b1Map[b1Map < 0] <- 0
        if (any(dim(b1Map) != mpmESTATICSModel$sdim)) 
            stop("dimension of B1 map does not match data dimension")
    }
    else {
        if (verbose) 
            cat("no B1 correction\n")
        b1Map <- array(1, mpmESTATICSModel$sdim)
    }
    t1FA <- mpmESTATICSModel$FA[1]
    pdFA <- mpmESTATICSModel$FA[length(mpmESTATICSModel$t1Files) + 
        length(mpmESTATICSModel$mtFiles) + 1]
    t1TR <- mpmESTATICSModel$TR[1]
    pdTR <- mpmESTATICSModel$TR[length(mpmESTATICSModel$t1Files) + 
        length(mpmESTATICSModel$mtFiles) + 1]
    if (verbose) 
        cat("calculating R1 ... ")
    alphat1 <- b1Map * t1FA/180 * pi
    alphapd <- b1Map * pdFA/180 * pi
    SINalphat1 <- sin(alphat1)
    COSalphat1 <- cos(alphat1)
    SINalphapd <- sin(alphapd)
    COSalphapd <- cos(alphapd)
    rm(alphat1, alphapd)
    if (mpmESTATICSModel$model == 2) {
        enum <- mpmESTATICSModel$modelCoeff[1, , , ] - SINalphat1/SINalphapd * 
            mpmESTATICSModel$modelCoeff[3, , , ]
        denom <- mpmESTATICSModel$modelCoeff[1, , , ] * COSalphat1 - 
            SINalphat1/SINalphapd * mpmESTATICSModel$modelCoeff[3, 
                , , ] * COSalphapd
    }
    else {
        enum <- mpmESTATICSModel$modelCoeff[1, , , ] - SINalphat1/SINalphapd * 
            mpmESTATICSModel$modelCoeff[2, , , ]
        denom <- mpmESTATICSModel$modelCoeff[1, , , ] * COSalphat1 - 
            SINalphat1/SINalphapd * mpmESTATICSModel$modelCoeff[2, 
                , , ] * COSalphapd
    }
    E1 <- enum/denom
    rm(enum, denom, COSalphapd, SINalphapd)
    R1 <- -log(E1)/t1TR
    P2_a = getPolynomsP2_ab(pdTR, t1TR, pdFA, t1FA, verbose)$P2_a
    P2_b = getPolynomsP2_ab(pdTR, t1TR, pdFA, t1FA, verbose)$P2_b
    R1 = R1/((P2_a[1] * b1Map^2 + P2_a[2] * b1Map + P2_a[3]) * 
        R1 + (P2_b[1] * b1Map^2 + P2_b[2] * b1Map + P2_b[3]))
    E1 = exp(-R1 * t1TR)
    if (verbose) 
        cat("done\n")
    if (verbose) 
        cat("calculating PD ... ")
    enum <- (1 - COSalphat1 * E1) * mpmESTATICSModel$modelCoeff[1, 
        , , ] * mpmESTATICSModel$dataScale
    denom <- SINalphat1 * (1 - E1)
    PD <- enum/denom
    rm(enum, denom, SINalphat1)
    if (verbose) 
        cat("done\n")
    if (mpmESTATICSModel$model == 2) {
        if (verbose) 
            cat("calculating MT ... ")
        mtFA <- mpmESTATICSModel$FA[length(mpmESTATICSModel$t1Files) + 
            1]
        mtTR <- mpmESTATICSModel$TR[length(mpmESTATICSModel$t1Files) + 
            1]
        alphamt <- b1Map * mtFA/180 * pi
        E1mt <- E1^(mtTR/t1TR)
        E2mt <- E1^(TR2/t1TR)
        enom <- mpmESTATICSModel$modelCoeff[2, , , ] * mpmESTATICSModel$dataScale - 
            (1 - E2mt) * sin(alphamt) * PD
        denom <- mpmESTATICSModel$modelCoeff[2, , , ] * mpmESTATICSModel$dataScale * 
            cos(alphamt) * E1mt + PD * (E2mt - E1mt) * sin(alphamt)
        delta <- 1 - enom/denom
        rm(alphamt, enom, denom)
        delta = 100 * delta * (1 - 0.4)/(1 - 0.4 * b1Map)/b1Map^2
        if (verbose) 
            cat("done\n")
    }
    else {
        delta <- NULL
    }
    invisible(list(b1Map = b1Map, R1 = R1 * 1000, R2star = if (mpmESTATICSModel$model == 
        2) 1000 * mpmESTATICSModel$modelCoeff[4, , , ]/mpmESTATICSModel$TEScale else 1000 * 
        mpmESTATICSModel$modelCoeff[3, , , ]/mpmESTATICSModel$TEScale, 
        PD = PD, delta = delta, model = mpmESTATICSModel$model, 
        t1Files = mpmESTATICSModel$t1Files, mtFiles = mpmESTATICSModel$mtFiles, 
        pdFiles = mpmESTATICSModel$pdFiles, mask = mpmESTATICSModel$mask))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
