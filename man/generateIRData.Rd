\name{generateIRData}
\alias{generateIRData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generate IR MRI example data
}
\description{
The function generates  IR MRI example data for specified parameters
}
\usage{
generateIRData(segm, pCSF, pGM, pWM, InvTimes, sigma = 40)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segm}{
     array containing segmentation resuts for an 3D MRI template. Contains 1 for CSF,
     2 for Gray Matter and 3 for White Matter
}
  \item{pCSF}{
    Parameters (S,R) for CSF
}
  \item{pGM}{
    Parameters (f,R,S) for Gray Matter
}
  \item{pWM}{
    Parameters (f,R,S) for White Matter
}
  \item{InvTimes}{
    Vector of Inversion Times, length nTimes
}
  \item{sigma}{
    Noise standard variation
}
}
\value{
   array with dimension c(nTimes,dim(segm)) 
}
\author{
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{
 \code{\link{estimateIRfluid}}, \code{\link{estimateIRsolid}}, \code{\link{estimateIR}},\code{\link{smoothIRSolid}}, \code{\link{readIRData}}
}
\note{
   used in examples for IR functions
}
\keyword{ utililities }
 \concept{ IRMRI }
