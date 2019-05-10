\name{writeESTATICS}
\alias{writeESTATICS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write maps of ESTATICS parameters in standardized form as NIfTI files.
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
Object of class 'ESTATICSModel' or 'sESTATICSModel' as returned from function
\code{\link{estimateESTATICS}} or \code{\link{smoothESTATICS}}.}
  \item{dir}{
  Directory name (or path) for output.
}
  \item{prefix}{
  Prefix for file names
}
  \item{verbose}{
  logical - provide information on progress
}
}
\details{
   \code{R2}, \code{ST1}, \code{SPD} and, if available, \code{SMT}-maps are written
   as compressed NIfTI files into directory the speecified directory.
   If \code{class(mpmESTATICSModel) == "sESTATICSModel"} and an smoothed data are stored in
   \code{mpmESTATICSModel$smoothedData} the smoothed data are stored as ompressed NIfTI files
   in \code{dir} with filenames assembled using \code{prefix} and the names of the data source
   files.
}
\value{
The function returns \code{NULL}
}
\author{
Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
J\"org Polzehl \email{polzehl@wias-berlin.de}
}

\seealso{
\code{\link{readMPMData}}, \code{\link{estimateESTATICS}},
\code{\link{smoothESTATICS}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }% use one of  RShowDoc("KEYWORDS")
