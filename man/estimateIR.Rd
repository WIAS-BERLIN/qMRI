\name{estimateIR}
\alias{estimateIR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate IRMRI parameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimateIR(IRdata, InvTimes, segments, fixed = TRUE, smoothMethod = c("Depth", "PAWS"), bw = 5, TEScale = 100, dataScale = 1000, method = c("NLR", "QL"), sigma = NULL, L = 1, maxR2star = 50, varest = c("RSS", "data"), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IRdata}{
%%     ~~Describe \code{IRdata} here~~
}
  \item{InvTimes}{
%%     ~~Describe \code{InvTimes} here~~
}
  \item{segments}{
%%     ~~Describe \code{segments} here~~
}
  \item{fixed}{
%%     ~~Describe \code{fixed} here~~
}
  \item{smoothMethod}{
%%     ~~Describe \code{smoothMethod} here~~
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
}
  \item{TEScale}{
%%     ~~Describe \code{TEScale} here~~
}
  \item{dataScale}{
%%     ~~Describe \code{dataScale} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{maxR2star}{
%%     ~~Describe \code{maxR2star} here~~
}
  \item{varest}{
%%     ~~Describe \code{varest} here~~
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
function (IRdata, InvTimes, segments, fixed = TRUE, smoothMethod = c("Depth", 
    "PAWS"), bw = 5, TEScale = 100, dataScale = 1000, method = c("NLR", 
    "QL"), sigma = NULL, L = 1, maxR2star = 50, varest = c("RSS", 
    "data"), verbose = TRUE) 
{
    ergsFluid <- estimateIRfluid(IRdata, InvTimes, segments)
    Sfluid <- median(ergsFluid$Sfluid)
    Rfluid <- median(ergsFluid$Rfluid)
    ergsBrain <- erstimateIRsolid(IRdata, InvTimes, segments, 
        Sfluid, Rfluid)
    if (fixed) {
        if (smmothMethod == "Depth") 
            ergsSmooth <- SdepthSmooth(ergsBrain, segments)
        if (smmothMethod == "PAWS") 
            ergsSmooth <- vpawsaws::vpawscov2(ergsBrain, segments)
    }
    ergsSmooth
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
