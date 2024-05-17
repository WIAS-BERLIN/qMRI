\name{estimateIR}
\alias{estimateIR}
\title{
Estimate IRMRI parameters
}
\description{
Parameter estimation (intensity, relaxation rate, proportion of fluid) in Inversion Recovery MRI data.
}
\usage{
estimateIR(IRdataobj, TEScale = 100, dataScale = 1000, method = c("NLR", "QL"),
           varest = c("RSS","data"), fixed = TRUE, smoothMethod=c("PAWS","Depth"),
           kstar = 24, alpha = .025, bysegment = TRUE, verbose = TRUE)
}
\arguments{
  \item{IRdataobj}{
4D array of IRMRI signals. First dimension corresponds to Inversion times (InvTime).
}
  \item{TEScale}{
 Internal scale factor for Echo Times. This influences parameter scales in numerical calculations.

}
  \item{dataScale}{
 Internal scale factor for MR signals. This influences parameter scales in numerical calculations.
}
  \item{method}{
  Either \code{"NLS"} for nonlinear least squares (ignores Rician bias) or \code{"QL"} for Quasi-Likelihood. The second option is more accurate but requires additional information and is computationally more expensive.
}
  \item{varest}{
  Method to, in case of \code{method="QR"}, estimate \code{sigma}if not provided. Either from residual sums of squares (\code{"RSS"}) or MR signals (\code{"data"}) using function \code{varest} from package \code{aws}. Only to be used in case that no image registration was needed as preprocessing.
}
\item{fixed}{
 Should adaptive smoothing performed for Sx and Rx maps and fx maps reestimated afterwards ?
}
  \item{smoothMethod}{
 Either "PAWS" or "Depth". the second option is not yet implemented.
}

  \item{kstar}{
  number of steps used in PAWS
}
  \item{alpha}{
   significance level for decisions in aws algorithm (suggestion: between 1e-5 and 0.025)
}
   \item{bysegment}{
   \code{TRUE}: restrict smoothing to segments from segmentation, \code{FALSE}: restrict smoothing to solid mask.
}
  \item{verbose}{
Logical. Provide some runtime diagnostics.
}
}
\details{
This function implements the complete pipeline of IRMRI anlysis.
}
\value{
List of class \code{"IRmixed"} with components
\item{IRdata }{4D array containing the IRMRI data, first dimension refers to inversion times}
\item{InvTimes }{vector of inversion times}
\item{segm}{segmentation codes, 1 for CSF, 2 for GM, 3 for WM, 0 for out of brain}
\item{sigma}{noise standard deviation, if not specified estimated fron CSF areas in image with largest inversion time}
\item{L}{effective number of coils}
\item{fx}{Array of fluid proportions}
\item{Sx}{Array of maximal signals}
\item{Rx}{Array of relaxation rates}
\item{Sf}{Global estimate of maximal fluid signal}
\item{Rf}{Global estimate of fluid relaxation rate}
\item{ICovx}{Covariance matrix of estimates \code{fx}, \code{Sx} and \code{Rx}.}
\item{sigma}{Array of provided or estimated noise standard deviations }
\item{Convx}{Array of convergence indicators}
\item{rsdx}{Residual standard deviations}
The arrays contain entries for all voxel with \code{segments\%in\%1:3}.
}
\references{
J. Polzehl and K. Tabelow (2023), Magnetic Resonance Brain Imaging:
Modeling and Data Analysis Using R, 2nd Edition, Chapter 7, Springer, Use R! Series.
    <doi:10.1007/978-3-031-38949-8_7>.

J. Polzehl and K. Tabelow (2023), Magnetic Resonance Brain Imaging - 
Modeling and Data Analysis Using R: Code and Data.
    <doi:10.20347/WIAS.DATA.6>.

}
\author{
Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
J\"org Polzehl \email{polzehl@wias-berlin.de}

}
\seealso{
 \code{\link{estimateIRfluid}}, \code{\link{estimateIRsolid}}, \code{\link{estimateIRsolidfixed}},\code{\link{smoothIRSolid}}
}
\examples{
\donttest{
## runs about 30 seconds
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
InvTimes0 <- c(100, 200, 400, 600, 800, 1200, 1600, 2000, 2500, 3000, 
              3500, 4000, 4500, 5000, 6000, 15000)
nTimes <- length(InvTimes0)
sigma <- 40
## generate IR signal
IRdata <- generateIRData(segm, c(Sf,Rf), c(fgm,Rgm,Sgm), c(fwm,Rwm,Swm), InvTimes0, sigma)
for(i in 1:9) writeNIfTI(as.nifti(IRdata[i,,,]), 
                         file.path(dataDir,paste0("IR0",i)))
for(i in 10:nTimes) writeNIfTI(as.nifti(IRdata[i,,,]), 
                         file.path(dataDir,paste0("IR",i)))
## generate IRdata object
t1Files <- list.files(dataDir,"*.nii.gz",full.names=TRUE)
segmFile <- file.path(dataDir0,"Brainweb_segm")
IRdata <- readIRData(t1Files, InvTimes0, segmFile, sigma=sigma,
                     L=1, segmCodes=c("CSF","GM","WM"))
## estimate all
sIRmix <- estimateIR(IRdata, method="QL")
}

}
 \keyword{ models }
 \keyword{ regression }
 \concept{ IRMRI }
