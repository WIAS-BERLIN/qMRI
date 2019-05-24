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
      double precision function lkern(kern,xsq)
      implicit none
      integer kern
      double precision xsq,z
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
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vaws2(y,mask,nv,n1,n2,n3,hakt,lambda,theta,s2,bi,
     1                thnew,s2new,ncores,lwght,wght,swjy)
C
C   y        observed values of regression function
C   nv       number of vector components
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   kritical value
C   theta    estimates from last step   (input)
C   si2      assumes same variance in components and independence
C   bi       \sum  Wi   (output)
C   thnew    \sum  Wi Y / bi     (output)
C   ncores   number of cores
C
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit none

      integer nv,n1,n2,n3,ncores
      logical aws, mask(*)
      double precision y(nv,*),theta(nv,*),bi(*),thnew(nv,*),s2(*),
     1  lambda,wght(2),hakt,lwght(*),swjy(nv,ncores),s2new(*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,k,thrednr
      double precision bii,biinv,sij,swj,swj2,z,z1,z2,z3,wj,hakt2,
     1       hmax2,w1,w2,spmb,spf
      external lkern,KLdist2
      double precision lkern,KLdist2
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
C$OMP& ih3,lwght,wght,y,swjy,mask,s2new,s2)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,bii,biinv,swj,spmb,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,z,thrednr,swj2)
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
         swj2=0.d0
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
                     sij=KLdist2(theta(1,iind),theta(1,jind),s2(iind))
                     IF (sij.ge.biinv) CYCLE
                     IF (sij.gt.spmb) wj=wj*(1.d0-spf*(bii*sij-0.25d0))
                  END IF
                  swj=swj+wj
                  swj2=swj2+wj+wj
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
         z1=y(1,iind)-thnew(1,iind)
         z2=y(2,iind)-thnew(2,iind)
         if(swj.gt.1.1d0) THEN
            s2new(iind) = (z1*z1+z2*z2)/2.d0*swj2/(swj2-1.d0)
         ELSE
            s2new(iind) = s2(iind)
         END IF
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thnew,bi)
      RETURN
      END
      double precision function KLdistsi(thi,thj,si2,nv)
      implicit none
      integer nv
      double precision thi(nv), thj(nv), si2(nv,nv)
      integer k,l
      double precision z,zd,thik,zdk
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

      double precision function KLdist2(thi,thj,s2)
      implicit none
      double precision thi(2), thj(2), s2, z1, z2
      z1 = thi(1)-thj(1)
      z2 = thi(2)-thj(2)
      KLdist2 = (z1*z1+z2*z2)/s2
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
      double precision x,y,value,wght(2),eps,bw
      double precision fw1,fw2,fw3,z
      double precision sofw
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
      double precision function sofw(bw,kern,wght)
      implicit logical(a-z)
      integer kern
      double precision bw,wght(2)
      integer j1,j2,j3,dlw1,dlw2,dlw3,clw1,clw2,clw3,ih1,ih2,ih3
      double precision sw,sw2,h2,lkern,z1,z2,z3,z
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
      subroutine paramw3(h,vext,ind,w,n)
C  compute a description of local weights
C  h    - bandwidth
C  vext - vector (length 2) of relative voxel extensions
C  ind  - integer array dim (3,n) containing relative indices in xyz
C  w    - vector of corresponding weights
C  n    - number of positive weights (initial value
C         (2 int(h)+1)*(2 int(h/vext(1))+1)*(2 int(h/vext(2))+1)
      integer n,ind(3,n)
      double precision h,vext(2),w(n)
      integer i,i1,i2,i3,ih1,ih2,ih3
      double precision hsq,z1,z2,z3
      hsq=h*h
      ih1 = int(h)
      ih2 = int(h/vext(1))
      ih3 = int(h/vext(2))
      i=1
      DO i1=-ih1,ih1
         z1=i1*i1
         DO i2=-ih2,ih2
            z2=i2*vext(1)
            z2=z1+z2*z2
            IF(z2.ge.hsq) CYCLE
            DO i3=-ih3,ih3
               z3=i3*vext(2)
               z3=z2+z3*z3
               IF(z3.ge.hsq) CYCLE
               ind(1,i)=i1
               ind(2,i)=i2
               ind(3,i)=i3
               w(i)=1.d0-z3/hsq
               i=i+1
            END DO
         END DO
      END DO
      n=i-1
      RETURN
      END
      subroutine mediansm(y,mask,n1,n2,n3,ind,nind,work,ncores,yout)
C
C
C   3D median smoother of y with neighborhood defined by ind
C   results in yout
C   size of work needs to be 2*nind
C
      implicit none
      integer n1,n2,n3,nind,ind(3,nind),ncores
      logical mask(n1,n2,n3)
      double precision y(n1,n2,n3),yout(n1,n2,n3),work(nind,ncores)
      integer i1,i2,i3,j1,j2,j3,j,k,thrednr
      double precision fmedian
      external fmedian
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      thrednr = 1
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i1,i2,i3,j1,j2,j3,j,k,thrednr)
C$OMP DO SCHEDULE(GUIDED)
      DO i1=1,n1
!$         thrednr = omp_get_thread_num()+1
         DO i2=1,n2
            DO i3=1,n3
               if(.not.mask(i1,i2,i3)) THEN
                  yout(i1,i2,i3) = y(i1,i2,i3)
                  CYCLE
               ENDIF
               k=0
               DO j=1,nind
                  j1=i1+ind(1,j)
                  if(j1.le.0.or.j1.gt.n1) CYCLE
                  j2=i2+ind(2,j)
                  if(j2.le.0.or.j2.gt.n2) CYCLE
                  j3=i3+ind(3,j)
                  if(j3.le.0.or.j3.gt.n3) CYCLE
                  if(.not.mask(j1,j2,j3)) CYCLE
                  if(y(j1,j2,j3).le.0.d0) CYCLE
                  k=k+1
                  work(k,thrednr)=y(j1,j2,j3)
               END DO
               IF(k.gt.1) THEN
                  yout(i1,i2,i3) = fmedian(work(1,thrednr),k)
C                  call qsort3(work(1,thrednr),1,k)
C                  IF (mod(k,2) == 0) THEN
C                     yout(i1,i2,i3) =
C     1               (work(k/2,thrednr)+work(k/2+1,thrednr))/2.d0
C                  ELSE
C                     yout(i1,i2,i3) = work(k/2+1,thrednr)
C                  END IF
               ELSE
                  yout(i1,i2,i3) = y(i1,i2,i3)
               END IF
            END DO
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(yout)
      return
      end
      double precision function fmedian(x,n)
C
C  compute the median using a select algorithm instead of sorting
C  used in mediansm
C
      implicit none
      integer n
      double precision x(n)
      integer m
      m = n/2+1
      call qselect(x,n,m)
      fmedian = x(m)
      if(mod(n,2).eq.0) THEN
         m = n-m+1
         call qselect(x,n,m)
         fmedian = (fmedian+x(m))/2.d0
      END IF
      return
      end
      subroutine qselect(x,n,k)
C
C  partial sorting using a select algorithm
C  output x(k) contains the kth element of sort(x)
C
      implicit none
      integer k,n
      double precision x(n)
      integer lft,rght,cur,i
      double precision guess
      lft=1
      rght=n
      DO WHILE (lft.lt.rght)
         guess = x(k)
         call swap(x,k,rght)
         cur = lft
         DO i=cur,rght-1
            IF(x(i).lt.guess) THEN
               call swap(x,i,cur)
               cur=cur+1
            END IF
         END DO
         call swap(x,rght,cur)
         if(cur.eq.k) EXIT
         if(cur.lt.k) THEN
            lft=cur+1
         ELSE
            rght=cur-1
         END IF
      END DO
      RETURN
      END
      subroutine swap(x,k,l)
      implicit none
      integer k,l
      double precision x(*),t
      t = x(k)
      x(k) = x(l)
      x(l) = t
      return
      end
