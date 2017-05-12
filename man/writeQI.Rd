\name{writeQI}
\alias{writeQI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write QI parameter maps of to Nifti.

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
writeQI(qi, dir = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qi}{
%%     ~~Describe \code{qi} here~~
}
  \item{dir}{
%%     ~~Describe \code{dir} here~~
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
function (qi, dir = NULL, verbose = TRUE) 
{
    if (!is.null(dir)) {
        dir.create(dir)
        r2file <- file.path(dir, "R2")
        r1file <- file.path(dir, "R1")
        pdfile <- file.path(dir, "PD")
        mtfile <- file.path(dir, "MT")
    }
    else {
        r2file <- "R2"
        r1file <- "R1"
        pdfile <- "PD"
        mtfile <- "MT"
    }
    ds <- readNIfTI(qi$t1Files[1], reorient = FALSE)
    ds@datatype <- 16
    ds@magic <- "n+1"
    ds@vox_offset <- 352
    if (verbose) 
        cat("writing R2 file ... ")
    ds@descrip <- "R2"
    writeNIfTI(as.nifti(qi$R2star, ds), file = r2file)
    if (verbose) 
        cat("done\n")
    if (verbose) 
        cat("writing R1 file ... ")
    ds@descrip <- "R1"
    writeNIfTI(as.nifti(qi$R1, ds), file = r1file)
    if (verbose) 
        cat("done\n")
    if (verbose) 
        cat("writing PD file ... ")
    ds@descrip <- "PD"
    writeNIfTI(as.nifti(qi$PD, ds), file = pdfile)
    if (verbose) 
        cat("done\n")
    if (qi$model == 2) {
        if (verbose) 
            cat("writing MT file ... ")
        ds@descrip <- "MT"
        writeNIfTI(as.nifti(qi$delta, ds), file = mtfile)
        if (verbose) 
            cat("done\n")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
