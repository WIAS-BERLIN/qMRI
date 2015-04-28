#####################################################################################
#
#  local variance estimation for voxel within a mask assuming homogeneity with the mask
#
#####################################################################################
homogensdev <- function(img3d,mask,bw=5){
dimg <- dim(img3d)
z <- .Fortran("homgnsd",
               as.double(img3d),
               as.logical(mask),
               as.integer(dimg[1]),
               as.integer(dimg[2]),
               as.integer(dimg[3]),
               as.integer(bw),
               sd=double(prod(dimg)),
               ni=integer(prod(dimg)),
               PACKAGE="dti")[c("sd","ni")]
list(sd=array(z$sd,dimg),ni=array(z$ni,dimg))
}
####################################################################################
#
#   False discovery rate (on a mask)
#
####################################################################################
fdr <- function(pval, alpha=0.05, mask){
  pval[!mask] <- 1
  n <- length(pval[mask])
  ind <- 1:length(pval)
  oind <- order(pval)
  nind <- length(ind[pval[oind] <= ind/n*alpha])
  detected <- array(FALSE,dim(pval))
  detected[oind[1:nind]] <- TRUE
  cat("fdr threshold:",pval[oind[nind]],"\n")
  detected
}
##################################################################################
#
#  reassign voxel from cortex to white matter based on DIR, T1 and T1postKM
#  and classification of surrounding voxel
#
##################################################################################
reassignCortex2WM <- function(ind, indre, n0=2){
##  indwm     - extended WM (or cortex) definition
##  indre     - voxel that are probably not in cortex (or WM)
##  n0        - number of neighboring voxel that need to be in ind (out of 6) 
   ddim <- dim(ind)
   indwm <- .Fortran("rawmsegm",
                     indwm=as.logical(ind),
                     as.logical(indre),
                     as.integer(ddim[1]),
                     as.integer(ddim[2]),
                     as.integer(ddim[3]),
                     as.integer(n0),
                     PACKAGE="qMRI")$indwm
   array(indwm,ddim)
   }
################################################################################
#
#  detect border regions
#
################################################################################
detectborder <- function(simg,dw=1){
   ddim <- dim(simg)
   border <- .Fortran("detctbrd",
                      as.integer(simg),
                      as.integer(ddim[1]),
                      as.integer(ddim[2]),
                      as.integer(ddim[3]),
                      as.integer(dw),
                      border=logical(prod(ddim)),
                      PACKAGE="qMRI")$border
   array(border,ddim)
}

