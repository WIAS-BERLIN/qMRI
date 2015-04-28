C
C    Copyright (C) 2015 Weierstrass-Institut fuer 
C                       Angewandte Analysis und Stochastik (WIAS)
C
C    Author:  Joerg Polzehl
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.
C
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
C  USA.
C
C  The following routines are part of the aws package and contain  
C  FORTRAN 77 code needed in R functions aws, vaws, 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    location penalty for multivariate non-gridded aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function lmkern(kern,dx,xi,xj,h2)
      implicit logical (a-z)
      external lkern
      integer kern,dx,i
      real*8 xi(dx),xj(dx),h2,z,zd,lkern
      z=0.d0
      do 1 i=1,dx
         zd=xi(i)-xj(i)
         z=z+zd*zd
         if(z.gt.h2) goto 2
1     continue
      lmkern=lkern(kern,z/h2)
      goto 999
2     lmkern=0.d0
999   return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Local constant aws on a grid      
C
C   this is a reimplementation of the original aws procedure
C
C   should be slightly slower for non-Gaussian models (see function kldist)
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute the Kullback-Leibler Distance
C
C          Model=1    Gaussian   
C          Model=2    Bernoulli   
C          Model=3    Poisson   
C          Model=4    Exponential   
C          Model=5    Variance
C          Model=6    Noncentral Chi (Gaussian approximation, 
C                     variance mean dependence is introduces via factor bii)
C
C     computing dlog(theta) and dlog(1.d0-theta) outside the AWS-loops 
C     will reduces computational costs at the price of readability
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function kldist(model,thi,thj)
      implicit logical (a-z)
      integer model
      real*8 thi,thj,z,tthi
      IF (model.eq.1) THEN
C        Gaussian
         z=thi-thj
         kldist=z*z
      ELSE IF (model.eq.2) THEN
C        Bernoulli
         kldist=0.d0
         tthi=(1.d0-thi)
         IF (thi.gt.1.d-10) kldist=kldist+thi*log(thi/thj)
         IF (tthi.gt.1.d-10) kldist=kldist+tthi*log(tthi/(1.d0-thj))
      ELSE IF (model.eq.3) THEN
C        Poisson
         kldist=0.d0
         IF (thi.gt.1.d-10) kldist=thi*log(thi/thj)-thi+thj
      ELSE IF (model.eq.4) THEN
C        Exponential
         kldist=thi/thj-1.d0-log(thi/thj)
      ELSE IF (model.eq.5) THEN
C        Variance
         kldist=thi/thj-1.d0-log(thi/thj)
      ELSE IF (model.eq.6) THEN
C        Noncentral Chi with Gaussian approximation
         z=thi-thj
         kldist=z*z
      ELSE
C        use Gaussian
         z=thi-thj
         kldist=z*z
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute Location Kernel (Compact support only, based on x^2
C                                   ignores scaling)
C
C          Kern=1     Uniform
C          Kern=2     Epanechnicov
C          Kern=3     Biweight
C          Kern=4     Triweight
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function lkern(kern,xsq)
      implicit logical (a-z)
      integer kern
      real*8 xsq,z
      IF (xsq.ge.1) THEN
         lkern=0.d0
      ELSE IF (kern.eq.1) THEN
         IF(xsq.le.0.5d0) THEN
            lkern=1.d0
         ELSE
            lkern=2.d0*(1.d0-xsq)
         END IF
      ELSE IF (kern.eq.2) THEN
         lkern=1.d0-xsq
      ELSE IF (kern.eq.3) THEN
         z=1.d0-xsq
         lkern=z*z
      ELSE IF (kern.eq.4) THEN
         z=1.d0-xsq
         lkern=z*z*z
      ELSE IF (kern.eq.5) THEN
         lkern=exp(-xsq*8.d0)
      ELSE
C        use Epanechnikov
         lkern=1.d0-xsq
      ENDIF
      RETURN 
      END   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        Compute truncated Exponential Kernel 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function skern(x,xmin,xmax)
      implicit logical (a-z)
      real*8 x,xmin,xmax,spf
      spf=xmax/(xmax-xmin)
      IF (x.le.xmin) THEN
         skern=1.d0
      ELSE IF (x.gt.xmax) THEN
         skern=0.d0
      ELSE
         skern=exp(-spf*(x-xmin))
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vaws(y,mask,nv,n1,n2,n3,hakt,lambda,theta,si2,bi,
     1                thnew,ncores,lwght,wght,swjy)
C
C   y        observed values of regression function 
C   nv       number of vector components
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   kritical value
C   theta    estimates from last step   (input)
C   si2      inverse covariance matrices of vectors in y 
C   bi       \sum  Wi   (output)
C   thnew    \sum  Wi Y / bi     (output)
C   ncores   number of cores
C   
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)

      integer nv,n1,n2,n3,ncores
      logical aws, mask(*)
      real*8 y(nv,*),theta(nv,*),bi(*),thnew(nv,*),si2(nv,nv,*),lambda,
     1  wght(2),hakt,lwght(*),swjy(nv,ncores)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,k,thrednr
      real*8 bii,biinv,sij,swj,z,z1,z2,z3,wj,hakt2,hmax2,w1,w2,spmb,spf
      external lkern,KLdistsi
      real*8 lkern,KLdistsi
!$      integer omp_get_thread_num 
!$      external omp_get_thread_num
      thrednr = 1
C just to prevent a compiler warning
      hakt2=hakt*hakt
      spf=4.d0/3.d0
      ih1=FLOOR(hakt)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      hmax2=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(2,(z1*z1+z2)/hakt2)
               if(lwght(jind).gt.0.d0) hmax2=max(hmax2,z2+z1*z1)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(thnew,bi,nv,n1,n2,n3,hakt2,hmax2,theta,
C$OMP& ih3,lwght,wght,y,swjy,mask)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& model,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,bii,biinv,swj,spmb,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,z,thrednr)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         if(.not.mask(iind)) CYCLE
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1         
C    nothing to do, final estimate is already fixed by control 
         bii=bi(iind)/lambda
         biinv=1.d0/bii
         spmb=0.25d0/bii
C   scaling of sij outside the loop
         swj=0.d0
         DO k=1,nv
            swjy(k,thrednr)=0.d0
         END DO
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  if(.not.mask(jind)) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  IF (aws) THEN
                     sij=KLdistsi(theta(1,iind),theta(1,jind),
     1                            si2(1,1,iind),nv)
                     IF (sij.ge.biinv) CYCLE
                     IF (sij.gt.spmb) wj=wj*(1.d0-spf*(bii*sij-0.25d0))
                  END IF
                  swj=swj+wj
                  DO k=1,nv
                     swjy(k,thrednr)=swjy(k,thrednr)+wj*y(k,jind)
                  END DO
               END DO
            END DO
         END DO
         DO k=1,nv
            thnew(k,iind)=swjy(k,thrednr)/swj
         END DO
         bi(iind)=swj
      END DO 
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thnew,bi)
      RETURN
      END
      real*8 function KLdistsi(thi,thj,si2,nv)
      implicit logical (a-z)
      integer nv
      real*8 thi(nv), thj(nv), si2(nv,nv)
      integer k,l
      real*8 z,zd,thik,zdk
      zd=thi(nv)-thj(nv)
      z=zd*zd*si2(nv,nv)
      DO k=1,nv-1
         thik=thi(k)
         zdk=thik-thj(k)
         z=z+zdk*zdk*si2(k,k)
         DO l=k+1,nv
            z=z+2.d0*(thi(l)-thj(l))*zdk*si2(k,l)
         END DO
      END DO
      KLdistsi=z
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vawsext(y,mask,nv,n1,n2,n3,yext,nve,hakt,lambda,
     1    theta,si2,bi,thnew,thext,ncores,lwght,wght,swjy,swjye)
C
C   y        observed values of regression function 
C   nv       number of vector components
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   kritical value
C   theta    estimates from last step   (input)
C   si2      inverse covariance matrices of vectors in y 
C   bi       \sum  Wi   (output)
C   thnew    \sum  Wi Y / bi     (output)
C   ncores   number of cores
C   
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)

      integer nv,n1,n2,n3,ncores,nve
      logical aws, mask(*)
      real*8 y(nv,*),theta(nv,*),bi(*),thnew(nv,*),si2(nv,nv,*),
     1       yext(nve,*),lambda,thext(nve,*),wght(2),hakt,lwght(*),
     2       swjy(nv,ncores),swjye(nve,ncores)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,k,thrednr
      real*8 bii,biinv,sij,swj,z,z1,z2,z3,wj,hakt2,hmax2,w1,w2,spmb,spf
      external lkern,KLdistsi
      real*8 lkern,KLdistsi
!$      integer omp_get_thread_num 
!$      external omp_get_thread_num
      thrednr = 1
C just to prevent a compiler warning
      hakt2=hakt*hakt
      spf=4.d0/3.d0
      ih1=FLOOR(hakt)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      hmax2=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(2,(z1*z1+z2)/hakt2)
               if(lwght(jind).gt.0.d0) hmax2=max(hmax2,z2+z1*z1)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(thnew,bi,nv,n1,n2,n3,hakt2,hmax2,theta,
C$OMP& ih3,lwght,wght,y,swjy,mask,swjye,yext,thext,nve)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& model,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,bii,biinv,swj,spmb,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,z,thrednr)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         if(.not.mask(iind)) CYCLE
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1         
C    nothing to do, final estimate is already fixed by control 
         bii=bi(iind)/lambda
         biinv=1.d0/bii
         spmb=0.25d0/bii
C   scaling of sij outside the loop
         swj=0.d0
         DO k=1,nv
            swjy(k,thrednr)=0.d0
         END DO
         DO k=1,nve
            swjye(k,thrednr)=0.d0
         END DO
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  if(.not.mask(jind)) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  IF (aws) THEN
                     sij=KLdistsi(theta(1,iind),theta(1,jind),
     1                            si2(1,1,iind),nv)
                     IF (sij.ge.biinv) CYCLE
                     IF (sij.gt.spmb) wj=wj*(1.d0-spf*(bii*sij-0.25d0))
                  END IF
                  swj=swj+wj
                  DO k=1,nv
                     swjy(k,thrednr)=swjy(k,thrednr)+wj*y(k,jind)
                  END DO
                  DO k=1,nve
                     swjye(k,thrednr)=swjye(k,thrednr)+wj*yext(k,jind)
                  END DO
               END DO
            END DO
         END DO
         DO k=1,nv
            thnew(k,iind)=swjy(k,thrednr)/swj
         END DO
         DO k=1,nve
            thext(k,iind)=swjye(k,thrednr)/swj
         END DO
         bi(iind)=swj
      END DO 
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thnew,bi,thext)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   determine sum of location weights for a given geometry a(3) and given 
C   bandwidth
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Algorithmus zur Nullstellenbestimmung einer monotonen Funktion auf(0,\infty)
      subroutine gethani(x,y,kern,value,wght,eps,bw)
      implicit logical(a-z)
      integer kern
      real*8 x,y,value,wght(2),eps,bw
      real*8 fw1,fw2,fw3,z
      real*8 sofw
      external sofw
      if(x.ge.y) RETURN
      fw1=sofw(x,kern,wght)
      fw2=sofw(y,kern,wght)
      DO WHILE(fw1.gt.value)
         x=x*x/y
         fw1=sofw(x,kern,wght)
      END DO
      DO WHILE(fw2.le.value)
         y=y*y/x
         fw2=sofw(y,kern,wght)
      END DO
      DO WHILE(min(fw2/value,value/fw1).gt.1.d0+eps)
         z=x+(value-fw1)/(fw2-fw1)*(y-x)
         fw3=sofw(z,kern,wght)
         if(fw3.le.value) THEN
            x=z
            fw1=fw3
         ENDIF
         if(fw3.ge.value) THEN
            y=z
            fw2=fw3
         ENDIF
         call rchkusr()
      END DO
      if(fw2/value.gt.value/fw1) THEN
          bw=x+(value-fw1)/(fw2-fw1)*(y-x)
      ELSE
          bw=y-(fw2-value)/(fw2-fw1)*(y-x)
      ENDIF
      RETURN
      END  
      real*8 function sofw(bw,kern,wght)
      implicit logical(a-z)
      integer kern
      real*8 bw,wght(2)
      integer j1,j2,j3,dlw1,dlw2,dlw3,clw1,clw2,clw3,ih1,ih2,ih3
      real*8 sw,sw2,h2,lkern,z1,z2,z3,z
      external lkern
      h2=bw*bw
C
C   first calculate location weights
C
      ih3=FLOOR(bw*wght(2))
      ih2=FLOOR(bw*wght(1))
      ih1=FLOOR(bw)
      dlw1=2*ih1+1
      dlw2=2*ih2+1
      dlw3=2*ih3+1
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
      sw=0.d0
      sw2=0.d0
      DO j1=1,dlw1
         z1=(clw1-j1)
         z1=z1*z1
         if(wght(1).gt.0.d0) THEN
            ih2=FLOOR(sqrt(h2-z1)*wght(1))
            DO j2=clw2-ih2,clw2+ih2
               z2=(clw2-j2)/wght(1)
               z2=z1+z2*z2
               if(wght(2).gt.0.d0) THEN
                  ih3=FLOOR(sqrt(h2-z2)*wght(2))
                  DO j3=clw3-ih3,clw3+ih3
                     z3=(clw3-j3)/wght(2)
                     z=lkern(kern,(z3*z3+z2)/h2)
                     sw=sw+z
                     sw2=sw2+z*z
                  END DO
               ELSE
                  z=lkern(kern,z2/h2)
                  sw=sw+z
                  sw2=sw2+z*z
               END IF
            END DO
         ELSE
            z=lkern(kern,z1/h2)
            sw=sw+z
            sw2=sw2+z*z
         END IF
      END DO
      sofw=sw*sw/sw2
      RETURN
      END
