rimage <- function(x = seq(0, 1, length.out = nrow(z)), 
                   y = seq(0, 1, length.out = ncol(z)),
                   z, zlim=NULL, col=grey(0:255/255),
                   low="blue", up="gold", NAcolor="red",...){
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        if (is.null(dim(x))) 
          stop("argument must be matrix-like")
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
        y <- seq.int(0, 1, length.out = ncol(z))
      }
    }
  }
   if(is.null(zlim)){ 
     zlim<-range(z)
   } else {
     z[z<zlim[1]] <- zlim[1]
     z[z>zlim[2]] <- zlim[2]
     col <- c(low,col,up)
   }
   image(x, y, z, zlim=zlim, col=col, ...)
   if(any(is.na(z))) image(x,y,is.na(z),col=c(NA,NAcolor),add=TRUE)
}
