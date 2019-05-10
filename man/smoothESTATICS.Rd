\name{smoothESTATICS}
\alias{smoothESTATICS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adaptive smoothing of ESTATICS parameters and MPM data
}
\description{
Performs adaptive smoothing of parameter maps in the ESTATICS model and if
\code{mpmData} is specified these data. Implements both vectorized variants of the
Adaptive Weights Smoothing (AWS, Polzehl and Spokoiny (2006))
and patchwise AWS (PAWS, Polzehl et al (2018)) algorithms with weighting schemes determined by
the estimated parameter maps and their covariances.
}
\usage{
smoothESTATICS(mpmESTATICSModel, mpmData = NULL, kstar = 16, alpha = 0.025,
               patchsize = 0, wghts = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mpmESTATICSModel}{
Object of class 'ESTATICSModel' as returned from function \code{\link{estimateESTATICS}}.
}
  \item{mpmData}{
(optional) Object of class MPMData as created by \code{\link{readMPMData}} from which
the parameter maps were obtained.
}
  \item{kstar}{
Maximum number of steps.
}
  \item{alpha}{
specifies the scale parameter for the adaptation criterion. smaller values are more restrictive.
}
  \item{patchsize}{
Patchsize in PAWS, 0 corresponds to AWS, alternative values are 1 and 2.
}
  \item{wghts}{
 (optional) voxel size if measurments are not isotropic.
}
  \item{verbose}{
  logical - provide information on progress
}
}
\value{
list with components
                 \item{modelCoeff }{Estimated parameter maps}
                 \item{invCov }{ map of inverse covariance matrices}
                 \item{isConv }{ convergence indicator map}
                 \item{bi}{ Sum of weights map from AWS/PAWS}
                 \item{smoothPar}{ smooting parameters used in AWS/PAWS}
                 \item{smoothedData}{ smoothed mpmData}
                 \item{sdim }{ image dimension}
                 \item{nFiles }{ number of images}
                 \item{t1Files }{ vector of T1 filenames}
                 \item{pdFiles }{ vector of PD filenames}
                 \item{mtFiles }{ vector of MT filenames}
                 \item{model }{ model used (depends on specification of MT files)}
                 \item{maskFile }{ filename of brain mask}
                 \item{mask }{ brain mask}
                 \item{sigma }{ sigma}
                 \item{L }{ L}
                 \item{TR }{ TR values}
                 \item{TE }{ TE values}
                 \item{FA }{ Flip angles (FA)}
                 \item{TEScale }{ TEScale}
                 \item{dataScale }{dataScale}
and class-attribute 'sESTATICSModel'
}
\references{
Joerg Polzehl, Vladimir Spokoiny,  Propagation-separation approach for local likelihood estimation,  Probab. Theory Related Fields  135 (3),  (2006)  , pp. 335--362.

Joerg Polzehl, Kostas Papafitsorus, Karsten Tabelow (2018). Patch-wise adaptive weights smoothing. WIAS-Preprint 2520.
}
\author{
Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
J\"org Polzehl \email{polzehl@wias-berlin.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{readMPMData}}, \code{\link{estimateESTATICS}}
}
\keyword{ model }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
