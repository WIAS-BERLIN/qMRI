\name{writeQI}
\alias{writeQI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write estimated maps  in standardized form as NIfTI files.
}
\description{
Quantitative \code{R2}, \code{R1}, \code{PD} and, if available, \code{MT}-maps are written
as compressed NIfTI files into directory the specified directory.

}
\usage{
writeQI(qi, dir = NULL, prefix="qmap", verbose = TRUE)
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
\item{prefix}{
Prefix for file names
}
  \item{verbose}{
logical - provide information on progress
}
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
\examples{
  # run with example{"writeQI",run.dontrun=TRUE}
\dontrun{
dataDir <- system.file("extdata",package="qMRI")
#
#  set file names for T1w, MTw and PDw images
#
t1Names <- paste0("t1w_",1:8,".nii.gz")
mtNames <- paste0("mtw_",1:6,".nii.gz")
pdNames <- paste0("pdw_",1:8,".nii.gz")
t1Files <- file.path(dataDir, t1Names)
mtFiles <- file.path(dataDir, mtNames)
pdFiles <- file.path(dataDir, pdNames)
#
#  file names of mask and B1 field map
#
B1File <- file.path(dataDir, "B1map.nii.gz")
maskFile <- file.path(dataDir, "mask.nii.gz")
#
#  Acquisition parameters (TE, TR, Flip Angle) for T1w, MTw and PDw images
#
TE <- c(2.3, 4.6, 6.9, 9.2, 11.5, 13.8, 16.1, 18.4,
        2.3, 4.6, 6.9, 9.2, 11.5, 13.8,
        2.3, 4.6, 6.9, 9.2, 11.5, 13.8, 16.1, 18.4)
TR <- rep(25, 22)
FA <- c(rep(21, 8), rep(6, 6), rep(6, 8))
#
#   read MPM example data
#
library(qMRI)
mpm <- readMPMData(t1Files, pdFiles, mtFiles,
                   maskFile, TR = TR, TE = TE,
                   FA = FA, verbose = FALSE)
#
#  Estimate Parameters in the ESTATICS model
#
modelMPM <- estimateESTATICS(mpm, method = "NLR")
# Alternatively using Quasi-Likelihood
sigma <- 50
modelMPMQL <- estimateESTATICS(mpm, method = "QL",
                  sigma = array(sigma,mpm$sdim), L = 1)
#
#  smooth maps of ESTATICS Parameters
#
setCores(2)
modelMPMQLsp1 <- smoothESTATICS(modelMPM,
                              kstar = 16,
                              alpha = 0.004,
                              patchsize=1,
                              verbose = TRUE)
#
#  resulting ESTATICS parameter maps for central coronal slice
#
if(require(adimpro)){
rimage.options(zquantiles=c(.01,.99), ylab="z")
par(mfrow=c(2,4),mar=c(3,3,3,1),mgp=c(2,1,0))
pnames <- c("T1","MT","PD","R2star")
for(i in 1:4){
   rimage(modelMPMQL$modelCoeff[i,,11,])
   title(pnames[i])
   }
for(i in 1:4){
   rimage(modelMPMQLsp1$modelCoeff[i,,11,])
   title(paste("smoothed",pnames[i]))
   }
}
#
#  Compute quantitative maps (R1, R2star, PD, MT)
#
qMRIQLMaps <- calculateQI(modelMPMQL,
                        b1File = B1File,
                        TR2 = 3.4)
qMRIQLSmoothedp1Maps <- calculateQI(modelMPMQLsp1,
                                    b1File = B1File,
                                    TR2 = 3.4)
#
#  resulting quantitative maps for central coronal slice
#
if(require(adimpro)){
rimage.options(zquantiles=c(.01,.99), ylab="z")
par(mfrow=c(2,4),mar=c(3,3,3,1),mgp=c(2,1,0))
nmaps <- c("R1","R2star","PD","MT")
qmap <- extract(qMRIQLMaps,nmaps)
for (i in 1:4) rimage(qmap[[i]][,11,],main=nmaps[i])
qmap <- extract(qMRIQLSmoothedp1Maps,nmaps)
for (i in 1:4) rimage(qmap[[i]][,11,],main=paste("Smoothed",nmaps[i]))
}
#
#  write qmaps
#
writeQI(qMRIQLSmoothedp1Maps,prefix="sqmap")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }% use one of  RShowDoc("KEYWORDS")
