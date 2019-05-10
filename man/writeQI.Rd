\name{writeQI}
\alias{writeQI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write estimated maps  in standardized form as NIfTI files.
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
  Object of class 'qMaps' as returned from function
  \code{\link{calculateQI}}
}
  \item{dir}{
Directory name (or path) for output.
}
  \item{verbose}{
logical - provide information on progress
}
}
\details{
Quantitative \code{R2}, \code{R1}, \code{PD} and, if available, \code{MT}-maps are written
as compressed NIfTI files into directory the speecified directory.
}
\value{
The function returns \code{NULL}
}
\author{
Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
J\"org Polzehl \email{polzehl@wias-berlin.de}
}

\seealso{
\code{\link{readMPMData}}, \code{\link{estimateESTATICS}},\code{\link{calculateQI}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }% use one of  RShowDoc("KEYWORDS")
