\name{MREdisplacement}
\alias{MREdisplacement}
\title{
Calculate the motion induced signal phase for IR-MRE in biphasic material 
}
\description{
The function takes magnitude images and phase images (as NIfTI files) recordet with inversion \code{IT1=Inf} 
and a second inversion time \code{IT2} that nulls the fluid signal. Tissue parameters (Relaxation rates)
are extracted from an object of class \code{"IRmixed"} calculated from data of a related IRMRI experiment. }
\usage{
MREdisplacement(MagnFiles1, PhaseFiles1, MagnFiles2, PhaseFiles2, TI2 = 2400,
                IRmixobj, method = c("full", "approx"),rescale=FALSE,verbose=FALSE)
}
\arguments{
  \item{MagnFiles1}{
  Filenames of magnitude images recorded with inversion time \code{IT=Inf} .
}
  \item{PhaseFiles1}{
  Filenames of phase images recorded with inversion time \code{IT=Inf} .
}
  \item{MagnFiles2}{
    Filenames of magnitude images recorded with inversion time \code{IT=IT2}.
}
  \item{PhaseFiles2}{
Filenames of phase images recorded with inversion time \code{IT=IT2} .
}
  \item{TI2}{
Inversion time used for \code{MagnFiles2} and \code{PhaseFiles2}. \code{IT2} should be selected to
extinguish the signal intendity for fluid.
}
  \item{IRmixobj}{
Object of class \code{"IRmixed"} obtained from a related IRMRI experiment.
}
  \item{method}{
Either \code{"full"} or \code{"approx"}
}
  \item{rescale}{
Logical, do we need to rescale phase images ? 
}
  \item{verbose}{
 Report scale range of phase images
}
}
\details{
The first 4 arguments need to be vectors of filenames of identical length with files containing compatible 3D NIfTI images. Object \code{IRmixobj} needs to contain a components \code{segm} and \code{Rx} of compatible dimension that need to be registered to the MRE images.
}
\value{
A list of class \code{"IRMREbiphasic"} with components
\item{phisolid}{displacement solid}
\item{phifluid}{displacement fluid}
}
\references{
J. Polzehl and K. Tabelow (2023), Magnetic Resonance Brain Imaging:
Modeling and Data Analysis Using R, 2nd Edition, Chapter 7, Springer, Use R! Series.
    <doi:10.1007/978-3-031-38949-8_7>.
}

\author{
Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
J\"org Polzehl \email{polzehl@wias-berlin.de}
}

\seealso{
 \code{\link{estimateIRfluid}}, \code{\link{estimateIRsolid}}, \code{\link{estimateIRsolidfixed}},\code{\link{smoothIRSolid}}
}
\concept{ IR-MRE }
