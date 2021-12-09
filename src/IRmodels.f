      subroutine IRfluid(th, invtime, n, fval, grad)
C
C  function values and gradients (2 parameters, S_f, R_f)
C
C  fval = th(1) * abs( 1 - 2*exp(-invtime*th(2))
C
      implicit logical (a-z)
      integer n
      double precision th(2),invtime(n),fval(n),grad(n,2)
      integer i
      double precision z1,z2,th1,th2
      th1=th(1)
      th2=th(2)
      DO i=1,n
         z1=exp(-invtime(i)*th2)
         z2=1.d0-2.d0*z1
         fval(i)=th1*abs(z2)
         grad(i,1)=abs(z2)
         grad(i,2)=-dsign(2.d0*z1*invtime(i),z2)
      END DO
      RETURN
      END

      subroutine IRmix(th, invtime, s0, t1, n, fval, grad)
C
C  function values and gradients (3 parameters, S_s, R_s, f)
C
C  fval = abs(s0*th(1)*(1-2*exp(-invtime*t1))+(1-th(1))*th(3)*(1-2*exp(-invtime*th2)))
C
      implicit logical (a-z)
      integer n
      double precision th(3),invtime(n),s0,t1,fval(n),grad(n,3)
      integer i
      double precision z1,z2,z3,th1,th2,th3,th13
      th1=th(1)
      th2=th(2)
      th3=th(3)
      th13=th3*(1.d0-th1)
      DO i=1,n
         z1=s0*(1.d0-2.d0*exp(-invtime(i)*t1))
         z2=exp(-invtime(i)*th2)
         z3=1.d0-2.d0*z2
         fval(i)=abs(s0*th1*z1+th13*z3)
         grad(i,1)=dsign((z1-th3*z3),th1*z1+th13*z3)
         grad(i,2)=-dsign(th13*2.d0*z2*invtime(i),th1*z1+th13*z3)
         grad(i,3)=dsign((1.d0-th1)*z3,th1*z1+th13*z3)
      END DO
      RETURN
      END

      subroutine IRmix0(th1, invtime, th2, th3, s0, t1, n, fval, grad)
C
C  function values and gradients (1 parameters, f)
C
C  fval = abs(s0*th(1)*(1-2*exp(-invtime*t1))+(1-th(1))*th(3)*(1-2*exp(-invtime*th(2))))
C
      implicit logical (a-z)
      integer n
      double precision th1,invtime(n),s0,t1,fval(n),grad(n),th2,th3
      integer i
      double precision z1,z2,z3,th13
      th13=th3*(1.d0-th1)
      DO i=1,n
         z1=1.d0-2.d0*exp(-invtime(i)*t1)
         z2=exp(-invtime(i)*th2)
         z3=1.d0-2.d0*z2
         fval(i)=s0*abs(th1*z1+(1.d0-th1)*th3*z3)
         grad(i)=dsign(s0*(z1-th3*z3),th1*z1*th13*z3)
      END DO
      RETURN
      END


