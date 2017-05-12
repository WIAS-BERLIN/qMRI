\name{writeESTATICS}
\alias{writeESTATICS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write parameter maps of the ESTATICS model to Nifti.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
writeESTATICS(mpmESTATICSModel, dir = NULL, prefix = "sm", verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mpmESTATICSModel}{
%%     ~~Describe \code{mpmESTATICSModel} here~~
}
  \item{dir}{
%%     ~~Describe \code{dir} here~~
}
  \item{prefix}{
%%     ~~Describe \code{prefix} here~~
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
function (mpmESTATICSModel, dir = NULL, prefix = "sm", verbose = TRUE) 
{
    if (!is.null(dir)) {
        dir.create(dir)
        r2file <- file.path(dir, "R2")
        st1file <- file.path(dir, "ST1")
        spdfile <- file.path(dir, "SPD")
        smtfile <- file.path(dir, "SMT")
    }
    else {
        r2file <- "R2"
        st1file <- "ST1"
        spdfile <- "SPD"
        smtfile <- "SMT"
    }
    ds <- readNIfTI(mpmESTATICSModel$t1Files[1], reorient = FALSE)
    ds@datatype <- 16
    ds@magic <- "n+1"
    ds@vox_offset <- 352
    if (mpmESTATICSModel$model == 2) {
        if (verbose) 
            cat("writing R2 file ... ")
        ds@descrip <- "R2"
        writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[4, , 
            , ], ds), file = r2file)
        if (verbose) 
            cat("done\n")
        if (verbose) 
            cat("writing ST1 file ... ")
        ds@descrip <- "ST1"
        writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[1, , 
            , ], ds), file = st1file)
        if (verbose) 
            cat("done\n")
        if (verbose) 
            cat("writing SPD file ... ")
        ds@descrip <- "SPD"
        writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[3, , 
            , ], ds), file = spdfile)
        if (verbose) 
            cat("done\n")
        if (verbose) 
            cat("writing SMT file ... ")
        ds@descrip <- "SMT"
        writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[2, , 
            , ], ds), file = smtfile)
        if (verbose) 
            cat("done\n")
    }
    else {
        if (verbose) 
            cat("writing R2 file ... ")
        ds@descrip <- "R2"
        writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[3, , 
            , ], ds), file = r2file)
        if (verbose) 
            cat("done\n")
        if (verbose) 
            cat("writing ST1 file ... ")
        ds@descrip <- "ST1"
        writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[1, , 
            , ], ds), file = st1file)
        if (verbose) 
            cat("done\n")
        if (verbose) 
            cat("writing SPD file ... ")
        ds@descrip <- "SPD"
        writeNIfTI(as.nifti(mpmESTATICSModel$modelCoeff[2, , 
            , ], ds), file = spdfile)
        if (verbose) 
            cat("done\n")
    }
    if (!is.null(mpmESTATICSModel$smoothedData)) {
        ii <- 1
        for (i in 1:length(t1Files)) {
            ds <- readNIfTI(t1Files[i], reorient = FALSE)
            ds@datatype <- 16
            ds@magic <- "n+1"
            ds@vox_offset <- 352
            if (!is.null(dir)) {
                fname <- file.path(dir, paste(prefix, basename(t1Files[i]), 
                  sep = ""))
            }
            else {
                fname <- paste(prefix, basename(t1Files[i]), 
                  sep = "")
            }
            if (verbose) 
                cat("writing", fname, "... ")
            writeNIfTI(as.nifti(mpmESTATICSModel$smoothedData[ii, 
                , , ], ds), file = fname)
            if (verbose) 
                cat("done\n")
            ii <- ii + 1
        }
        if (mpmESTATICSModel$model == 2) {
            for (i in 1:length(mtFiles)) {
                ds <- readNIfTI(mtFiles[i], reorient = FALSE)
                ds@datatype <- 16
                ds@magic <- "n+1"
                ds@vox_offset <- 352
                if (!is.null(dir)) {
                  fname <- file.path(dir, paste(prefix, basename(mtFiles[i]), 
                    sep = ""))
                }
                else {
                  fname <- paste(prefix, basename(mtFiles[i]), 
                    sep = "")
                }
                if (verbose) 
                  cat("writing", fname, "... ")
                writeNIfTI(as.nifti(mpmESTATICSModel$smoothedData[ii, 
                  , , ], ds), file = fname)
                if (verbose) 
                  cat("done\n")
                ii <- ii + 1
            }
        }
        for (i in 1:length(pdFiles)) {
            ds <- readNIfTI(pdFiles[i], reorient = FALSE)
            ds@datatype <- 16
            ds@magic <- "n+1"
            ds@vox_offset <- 352
            if (!is.null(dir)) {
                fname <- file.path(dir, paste(prefix, basename(pdFiles[i]), 
                  sep = ""))
            }
            else {
                fname <- paste(prefix, basename(pdFiles[i]), 
                  sep = "")
            }
            if (verbose) 
                cat("writing", fname, "... ")
            writeNIfTI(as.nifti(mpmESTATICSModel$smoothedData[ii, 
                , , ], ds), file = fname)
            if (verbose) 
                cat("done\n")
            ii <- ii + 1
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
