align3D <- function(MRI,sourcemat,targetmat,sourcedim,targetdim){
  ##
  ##  realign MRI from source-space to target-space
  ##
  R1 <- R2 <- matrix(0,3,3)
  R1 <- solve(sourcemat[,1:3])
  R2 <- targetmat[,1:3]
  off1 <- sourcemat[,4]
  off2 <- targetmat[,4]
  dims <- dim(MRI)
  if(any(sourcedim!=dims)) stop("incompatible dimensions")
  dimt <- targetdim
  array(.Fortran("align3d",
                 as.double(MRI),
                 as.double(R1),
                 as.double(R2),
                 as.double(off2-off1),
                 as.integer(dims[1]),
                 as.integer(dims[2]),
                 as.integer(dims[3]),
                 as.integer(dimt[1]),
                 as.integer(dimt[2]),
                 as.integer(dimt[3]),
                 Mt=double(prod(dimt)),
                 PACKAGE="qMRI")$Mt,dimt)
}
