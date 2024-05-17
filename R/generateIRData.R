generateIRData <- function(segm, pCSF, pGM, pWM, InvTimes, sigma=40){
Sf <- pCSF[1]
Rf <- pCSF[1]
fgm <- pGM[1]
Rgm <- pGM[2]
Sgm <- pGM[3]
fwm <- pWM[1]
Rwm <- pWM[2]
Swm <- pWM[3]
fintCSF <- IRhomogen(c(Sf,Rf),InvTimes)
fintGM <- IRmix2(c(fgm,Rgm,Sgm),InvTimes,Sf,Rf)
fintWM <- IRmix2(c(fwm,Rwm,Swm),InvTimes,Sf,Rf)
## generate data
nTimes <- length(InvTimes)
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
IRdata
}
