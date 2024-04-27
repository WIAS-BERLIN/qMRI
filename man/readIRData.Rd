\name{readIRData}
\alias{readIRData}
\title{
   Prepare IRMRI dataset
}
\description{
The function reads IRMRI images given as NIfTI files in t1Files, inversion times and segmentation image(s) aund prepares an object class \code{"IRdata"}
}
\usage{
readIRData(t1Files, InvTimes, segmFile, sigma = NULL, L = 1,
           segmCodes = c("GM", "WM", "CSF"))
}
\arguments{
  \item{t1Files}{
   Names of NIfTI files containing the recorded images.
}
  \item{InvTimes}{
   Corresponding inversion times
}
  \item{segmFile}{
   Either a NIfTI file containing a segmentation into GM, WM and CSF or
   three files containing probability maps for GM, WM and CSF
}
  \item{sigma}{
   Noise standard deviation
}
  \item{L}{
   Effective number of coils, L=1 assumes a Rician signal distribution
}
  \item{segmCodes}{
   sequence of tissue code in segmFile
}
}
\value{
A list of class \code{"IRdata"} with components
\item{IRdata }{4D array containing the IRMRI data, first dimension refers to inversion times}
\item{InvTimes }{vector of inversion times}
\item{segm}{segmentation codes, 1 for CSF, 2 for GM, 3 for WM, 0 for out of brain}
\item{sigma}{noise standard deviation, if not specified estimated fron CSF areas in image with largest inversion time}
\item{L}{effective number of coils}
}
\references{
J. Polzehl and K. Tabelow (2023), Magnetic Resonance Brain Imaging:
Modeling and Data Analysis Using R, 2nd Edition, Chapter 7, Springer, Use R! Series.
    <doi:10.1007/978-3-031-38949-8_7>.

J. Polzehl and K. Tabelow (2023), Magnetic Resonance Brain Imaging - 
Modeling and Data Analysis Using R: Code and Data.
    \url{https://doi.org/10.20347/WIAS.DATA.6}.

}
\author{
Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{
 \code{\link{estimateIRfluid}}, \code{\link{estimateIRsolid}}, \code{\link{estimateIR}},\code{\link{smoothIRSolid}}
}
\keyword{ utililities }
 \concept{ IRMRI }
