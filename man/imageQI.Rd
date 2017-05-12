\name{imageQI}
\alias{imageQI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Display a slice of 3D QI parameter maps.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
imageQI(qi, view = 1, slice = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qi}{
%%     ~~Describe \code{qi} here~~
}
  \item{view}{
%%     ~~Describe \code{view} here~~
}
  \item{slice}{
%%     ~~Describe \code{slice} here~~
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
function (qi, view = 1, slice = 1) 
{
    mask <- switch(view, qi$mask[slice, , ], qi$mask[, slice, 
        ], qi$mask[, , slice])
    r2star <- switch(view, qi$R2star[slice, , ], qi$R2star[, 
        slice, ], qi$R2star[, , slice])
    r2star[!mask] <- 0
    r1 <- switch(view, qi$R1[slice, , ], qi$R1[, slice, ], qi$R1[, 
        , slice])
    r1[!mask] <- 0
    pd <- switch(view, qi$PD[slice, , ], qi$PD[, slice, ], qi$PD[, 
        , slice])
    pd[!mask] <- 0
    if (qi$model == 2) {
        delta <- switch(view, qi$delta[slice, , ], qi$delta[, 
            slice, ], qi$delta[, , slice])
        delta[!mask] <- 0
    }
    indx <- 1:dim(r2star)[1]
    indy <- 1:dim(r2star)[2]
    if (qi$model == 2) {
        def.par <- par(mfrow = c(2, 2), mar = c(3, 3, 3, 0))
        rimage(indx, indy, r2star, zlim = c(0, 0.05), main = "R2star")
        rimage(indx, indy, r1, zlim = c(2e-04, 0.0015), main = "R1")
        rimage(indx, indy, pd, zlim = c(0, 10000), main = "PD")
        rimage(indx, indy, delta, zlim = c(0, 0.03), main = "MT")
    }
    else {
        def.par <- par(mfrow = c(2, 2), mar = c(3, 3, 3, 0))
        rimage(indx, indy, r2star, zlim = c(0, 0.05), main = "R2star")
        rimage(indx, indy, r1, zlim = c(2e-04, 0.0015), main = "R1")
        rimage(indx, indy, pd, zlim = c(0, 10000), main = "PD")
    }
    par(def.par)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
