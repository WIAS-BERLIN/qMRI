\name{readMPMData}
\alias{readMPMData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read datasets from quantitative MR imaging}

\description{
The function reads data generated in Multimodal Parameter Mapping (MPM) experiments.
}
\usage{
readMPMData(t1Files = NULL, pdFiles = NULL, mtFiles = NULL, maskFile = NULL, sdim = NULL, TR = NULL, TE = NULL, FA = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t1Files}{
    Vector of filenames corresponding to T1 weighted images (in Nifti-Format) with varying TE 
}
  \item{pdFiles}{
    Vector of filenames corresponding to PD weighted images (in Nifti-Format) with varying TE 
}
  \item{mtFiles}{
   optional Vector of filenames corresponding to MT weighted images (in Nifti-Format) with varying TE 
}
  \item{maskFile}{
   optional filename for mask
}
  \item{sdim}{
   array dimension of .nii files
}
  \item{TR}{
   optional TR vector, if omitted information is extracted from .nii files if possible
}
  \item{TE}{
   optional TE vector, if omitted information is extracted from .nii files if possible
}
  \item{FA}{
   optional FA (flip-angle) vector, if omitted information is extracted from .nii files if possible
}
  \item{verbose}{
   provide information on progress
}
}
\details{
}
\value{
List with components
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
     \item{ddata}{mpm data}
                 \item{sdim }{ sdim}
                 \item{nFiles }{ nFiles}
                \item{ t1Files }{ t1Files}
                \item{ pdFiles }{ pdFiles}
                \item{mtFiles }{ mtFiles}
                \item{ model }{ model}
                 \item{maskFile }{ maskFile}
                 \item{mask }{ mask}
                 \item{TR }{ TR}
                 \item{TE }{ TE}
                 \item{FA }{ FA}
}
\references{
Weiskopf, N.; Suckling, J.; Williams, G.; Correia, M. M.; Inkster, B.; Tait, R.; Ooi, C.; Bullmore, E. T. & Lutti, A. Quantitative multi-parameter mapping of R1, PD(*), MT, and R2(*) at 3T: a multi-center validation. Front Neurosci, Wellcome Trust Centre for Neuroimaging, UCL Institute of Neurology, University College London London, UK., 2013, 7, 95
}
\author{
Karsten Tabelow \email{karsten.tabelow@wias-berlin.de}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{estimateESTATICS}}, \code{\link{estimateESTATICSQL}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
