      subroutine qflashm1(th,des,n,fval,grad)
C
C  function values and gradients (4 parameters)
C
      implicit logical (a-z)
      integer n
      real*8 th(4),des(n,4),fval(n),grad(n,4)
      integer i
      real*8 z1,z2,z3,z4,sa,ca,fv,g1,g2,g3,d1,d2,d3,d4,th1,th2,th3,th4
      th1=th(1)
      th2=th(2)
      th3=th(3)
      th4=th(4)
      DO i=1,n
         d1=des(i,1)
         d2=des(i,2)
         d3=des(i,3)
         d4=des(i,4)
         z1=exp(-th2*d1)
         z2=exp(-th3*d3)
         sa=sin(d2)
         ca=cos(d2)
         z3=1.d0-z2
         z4=1.d0-ca*z2
         fv=th1*z1*sa*z3/z4
         g1=z1*sa*z3/z4
         g2=-d1*fv
         g3=th1*z1*sa*((1.d0-z3/z4*ca)*d3*z2/z4)
         if(d4.gt.0) THEN
            g1=g1*th4
            g2=g2*th4
            g3=g3*th4
            grad(i,4)=fv
            fv=fv*th4
         ELSE
            grad(i,4)=0.d0
         END IF
         fval(i)=fv
         grad(i,1)=g1
         grad(i,2)=g2
         grad(i,3)=g3
      END DO
      RETURN
      END
      subroutine qflashm0(th,des,n,fval)
C
C  function values only
C
      implicit logical (a-z)
      integer n
      real*8 th(3),des(3,n),fval(n)
      integer i
      real*8 z1,z2
      DO i=1,n
         z1=exp(-th(2)*des(1,i))
         z2=exp(-th(3)*des(3,i))
         fval(i)=th(1)*z1*dsin(des(2,i))*(1.d0-z2)/
     1                       (1.d0-dcos(des(2,i))*z2)
      END DO
      RETURN
      END
      subroutine qflashpl(th,des,n,fval,grad)
C
C  function values and gradients (4 parameters)
C
      implicit logical (a-z)
      integer n
      real*8 th(4),des(n,4),fval(n),grad(n,4)
      integer i
      real*8 z4,fv
      DO i=1,n
         z4=exp(-th(4)*des(i,4))
         if(des(i,1).gt.0) THEN
            fv=z4*th(1)
            grad(i,1)=z4
            grad(i,2)=0.d0
            grad(i,3)=0.d0
         END IF
         if(des(i,2).gt.0) THEN
            fv=z4*th(2)
            grad(i,1)=0.d0
            grad(i,2)=z4
            grad(i,3)=0.d0
         END IF
         if(des(i,3).gt.0) THEN
            fv=z4*th(3)
            grad(i,1)=0.d0
            grad(i,2)=0.d0
            grad(i,3)=z4
         END IF
         grad(i,4)=-des(i,4)*fv
         fval(i)=fv
      END DO
      RETURN
      END

      subroutine qflashp0(th,r2star,des,n,fval,grad)
C
C  function values and gradients (3 parameters)
C
      implicit logical (a-z)
      integer n
      real*8 th(3),r2star,des(n,4),fval(n),grad(n,3)
      integer i
      real*8 z4,fv
      DO i=1,n
         z4=exp(-r2star*des(i,4))
         if(des(i,1).gt.0) THEN
            fv=z4*th(1)
            grad(i,1)=z4
            grad(i,2)=0.d0
            grad(i,3)=0.d0
         END IF
         if(des(i,2).gt.0) THEN
            fv=z4*th(2)
            grad(i,1)=0.d0
            grad(i,2)=z4
            grad(i,3)=0.d0
         END IF
         if(des(i,3).gt.0) THEN
            fv=z4*th(3)
            grad(i,1)=0.d0
            grad(i,2)=0.d0
            grad(i,3)=z4
         END IF
         fval(i)=fv
      END DO
      RETURN
      END

      subroutine qflashpl2(th,des,n,fval,grad)
C
C  function values and gradients (3parameters)
C
      implicit logical (a-z)
      integer n
      real*8 th(3),des(n,3),fval(n),grad(n,3)
      integer i
      real*8 z3,fv
      DO i=1,n
         z3=exp(-th(3)*des(i,3))
         if(des(i,1).gt.0) THEN
            fv=z3*th(1)
            grad(i,1)=z3
            grad(i,2)=0.d0
         END IF
         if(des(i,2).gt.0) THEN
            fv=z3*th(2)
            grad(i,1)=0.d0
            grad(i,2)=z3
         END IF
         grad(i,3)=-des(i,3)*fv
         fval(i)=fv
      END DO
      RETURN
      END

      subroutine qflashpl3(th,des,n,fval,grad)
C
C  function values and gradients (2parameters)
C
      implicit logical (a-z)
      integer n
      real*8 th(2),des(n,2),fval(n),grad(n,2)
      integer i
      real*8 z2,fv
      DO i=1,n
         z2=exp(-th(2)*des(i,2))
         if(des(i,1).gt.0) THEN
            fv=z2*th(1)
            grad(i,1)=z2
         END IF
         grad(i,2)=-des(i,2)*fv
         fval(i)=fv
      END DO
      RETURN
      END

