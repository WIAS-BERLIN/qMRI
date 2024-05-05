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
\examples{
\dontrun{
dataDir0 <- system.file("extdataIR", package = "qMRI")
dataDir <- tempdir("IRdata")
library(oro.nifti)
library(qMRI)
segm <- readNIfTI(file.path(dataDir0,"Brainweb_segm"))
Sf <- 900
Rf <- 0.000285
Sgm <- 400
Rgm <- 0.00075
fgm <- .15
Swm <- 370
Rwm <- 0.0011
fwm <- .05
InvTimes <- c(100, 200, 400, 600, 800, 1200, 1600, 2000, 2500, 3000, 
              3500, 4000, 4500, 5000, 6000, Inf)
InvTimes0 <- c(100, 200, 400, 600, 800, 1200, 1600, 2000, 2500, 3000, 
              3500, 4000, 4500, 5000, 6000, 15000)
x <- seq(100,12000,10)
## generate IR signal
fintCSF <- qMRI:::IRhomogen(c(Sf,Rf),InvTimes0)
fintGM <- qMRI:::IRmix2(c(fgm,Rgm,Sgm),InvTimes0,Sf,Rf)
fintWM <- qMRI:::IRmix2(c(fwm,Rwm,Swm),InvTimes0,Sf,Rf)
## generate data
sigma <- 40
nTimes <- length(InvTimes0)
nCSF <- sum(segm==1)
nGM <- sum(segm==2)
nWM <- sum(segm==3)
IRdata <- array(0,c(nTimes,prod(dim(segm))))
IRdata[,segm==1] <- sqrt(rnorm(nTimes*nCSF,fintCSF,sigma)^2+
                         rnorm(nTimes*nCSF,0,sigma)^2)
IRdata[,segm==2] <- sqrt(rnorm(nTimes*nGM,fintGM,sigma)^2+
                         rnorm(nTimes*nGM,0,sigma)^2)
IRdata[,segm==3] <- sqrt(rnorm(nTimes*nWM,fintWM,sigma)^2+
                         rnorm(nTimes*nWM,0,sigma)^2)
dim(IRdata) <- c(nTimes,dim(segm))
for(i in 1:9) writeNIfTI(as.nifti(IRdata[i,,,]), 
                         file.path(dataDir,paste0("IR0",i)))
for(i in 10:nTimes) writeNIfTI(as.nifti(IRdata[i,,,]), 
                         file.path(dataDir,paste0("IR",i)))
## generate IRdata object
t1Files <- list.files(dataDir,"*.nii.gz",full.names=TRUE)
segmFile <- file.path(dataDir0,"Brainweb_segm")
IRdata <- readIRData(t1Files, InvTimes0, segmFile, sigma=sigma,
                     L=1, segmCodes=c("CSF","GM","WM"))
}
}
\keyword{ utililities }
 \concept{ IRMRI }
