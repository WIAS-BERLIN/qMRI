\name{calculateQI}
\alias{calculateQI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain quantitative maps from estimated ESTATICS parameters.
}
\description{
Quantitaive imaging parameters are calculated from the estimated parameters
in the ESTATICS model. This involves a correction for magnetic field inhomogeneities
if the information is provided in argument \code{b1File} and use of a second
of a second recovery delay TR2 in case of Dual-Exitation FLASH measurements
(Helms 2008).
}
\usage{
calculateQI(mpmESTATICSModel, b1File = NULL, TR2 = 0, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mpmESTATICSModel}{
  Object of class 'ESTATICSModel' as returned from function \code{\link{estimateESTATICS}}.
}
  \item{b1File}{
  (optional) Name of a file containing a B1-field inhomogeneity map (.nii)
}
  \item{TR2}{
   second recovery delay TR2 in case of Dual-Exitation FLASH measurements.
}
  \item{verbose}{
  logical: Monitor process.
}
}
\value{
List with components
        \item{b1Map }{ b1Map }
        \item{R1 }{ Estimated map of R1}
        \item{R2star }{ Estimated map of R2star}
        \item{PD }{Estimated map of  PD}
        \item{MT }{Estimated map of  delta (if MT-series was used)}
        \item{model }{Type of ESTATICS model used}
        \item{t1Files }{filenames T1}
        \item{mtFiles }{filenames MT}
        \item{pdFiles }{filenames PD}
        \item{mask }{brainmask}
and class-attribute 'qMaps'
.}
\references{
Helms, G.; Dathe, H.; Kallenberg, K. & Dechent, P.
High-Resolution Maps of Magnetization Transfer with Inherent Correction
for RF Inhomogeneity and T1 Relaxation Obtained from 3D FLASH MRI Magn.
Res. Med., 2008, 60, 1396-1407

Weiskopf, N.; Suckling, J.; Williams, G.; Correia, M. M.; Inkster, B.;
Tait, R.; Ooi, C.; Bullmore, E. T. & Lutti, A.
Quantitative multi-parameter mapping of R1, PD(*), MT, and R2(*) at 3T:
a multi-center validation. Front Neurosci,
Wellcome Trust Centre for Neuroimaging, UCL Institute of Neurology,
University College London, UK., 2013, 7, 95
}
\author{
Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
J\"org Polzehl \email{polzehl@wias-berlin.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{readMPMData}}, \code{\link{estimateESTATICS}},
\code{\link{smoothESTATICS}}, \code{\link{writeESTATICS}},
\code{\link[dti]{awslsigmc}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
