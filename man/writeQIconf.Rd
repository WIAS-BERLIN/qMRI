\name{writeQIconf}
\alias{writeQIconf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write confidence intervals of QI parameter maps of to Nifti.

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
writeQIconf(qiConf, dir = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qiConf}{
%%     ~~Describe \code{qiConf} here~~
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
function (qiConf, dir = NULL, verbose = TRUE) 
{
    if (!is.null(dir)) {
        dir.create(dir)
        r1file <- file.path(dir, "R1")
        r1Lfile <- file.path(dir, "R1lower")
        r1Ufile <- file.path(dir, "R1upper")
        r2file <- file.path(dir, "R2")
        r2Lfile <- file.path(dir, "R2lower")
        r2Ufile <- file.path(dir, "R2upper")
    }
    else {
        r1file <- "R1"
        r1Lfile <- "R1lower"
        r1Ufile <- "R1upper"
        r2file <- "R2"
        r2Lfile <- "R2lower"
        r2Ufile <- "R2upper"
    }
    ds <- readNIfTI(qiConf$t1Files[1], reorient = FALSE)
    ds@datatype <- 16
    ds@magic <- "n+1"
    ds@vox_offset <- 352
    if (verbose) 
        cat("writing R1 file ... ")
    ds@descrip <- "R1"
    writeNIfTI(as.nifti(qiConf$R1, ds), file = r1file)
    if (verbose) 
        cat("done\n")
    if (verbose) 
        cat("writing R1lower file ... ")
    ds@descrip <- "R1lower"
    writeNIfTI(as.nifti(qiConf$CIR1[1, , , ], ds), file = r1Lfile)
    if (verbose) 
        cat("done\n")
    if (verbose) 
        cat("writing R1upper file ... ")
    ds@descrip <- "R1upper"
    writeNIfTI(as.nifti(qiConf$CIR1[2, , , ], ds), file = r1Ufile)
    if (verbose) 
        cat("done\n")
    if (verbose) 
        cat("writing R2 file ... ")
    ds@descrip <- "R2"
    writeNIfTI(as.nifti(qiConf$R2star, ds), file = r2file)
    if (verbose) 
        cat("done\n")
    if (verbose) 
        cat("writing R2lower file ... ")
    ds@descrip <- "R2lower"
    writeNIfTI(as.nifti(qiConf$CIR2star[1, , , ], ds), file = r2Lfile)
    if (verbose) 
        cat("done\n")
    if (verbose) 
        cat("writing R2upper file ... ")
    ds@descrip <- "R2upper"
    writeNIfTI(as.nifti(qiConf$CIR2star[2, , , ], ds), file = r2Ufile)
    if (verbose) 
        cat("done\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
