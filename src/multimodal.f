      subroutine homgnsd(img,mask,n1,n2,n3,bw,sd,ni)
      implicit logical (a-z)
      integer n1,n2,n3,bw,ni(n1,n2,n3)
      real*8 img(n1,n2,n3),sd(n1,n2,n3)
      logical mask(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,anz
      real*8 z,z1,z2
C  calculate standard deviation for voxel within mask and with sup distance 
C of bw to the border
      DO i3=bw+1,n3-bw
         DO i2=bw+1,n2-bw
            DO i1=bw+1,n1-bw
               if(mask(i1,i2,i3)) THEN
                  z1=0.d0
                  z2=0.d0
                  anz=0
                  DO j3=i3-bw,i3+bw
                     DO j2=i2-bw,i2+bw
                        DO j1=i1-bw,i1+bw
                           if(mask(j1,j2,j3)) THEN
                              z=img(j1,j2,j3)
                              z1=z1+z
                              z2=z2+z*z
                              anz=anz+1
                           END IF
                        END DO
                     END DO
                  END DO
                  ni(i1,i2,i3)=anz-1
                  IF(anz.gt.1) THEN
                     z = z1/anz 
                     z = z2/anz - z*z
                     sd(i1,i2,i3) = sqrt(anz/(anz-1.d0)*z)
                  ELSE
                     sd(i1,i2,i3) = 0.d0
                  END IF
               ELSE
                  sd(i1,i2,i3) = 0.d0
                  ni(i1,i2,i3) = 0
               END IF
            END DO
         END DO
      END DO
      return
      end
C
C   reassign voxel between cortex and white matter
C
      subroutine rawmsegm(indwm,indre,n1,n2,n3,n0)
      implicit logical(a-z)
      integer n1,n2,n3,n0
      logical indwm(n1,n2,n3),indre(n1,n2,n3)
      integer i1,i2,i3,anz,anwm
      anz = 1
      DO WHILE (anz.gt.0)
         anz = 0
         DO i1=2,n1-1
            DO i2=2,n2-1
               DO i3=2,n3-1
                  IF(indre(i1,i2,i3)) THEN
                     anwm=0
C count WM neighbors
                     IF(indwm(i1-1,i2,i3)) anwm=anwm+1
                     IF(indwm(i1+1,i2,i3)) anwm=anwm+1
                     IF(indwm(i1,i2-1,i3)) anwm=anwm+1
                     IF(indwm(i1,i2+1,i3)) anwm=anwm+1
                     IF(indwm(i1,i2,i3-1)) anwm=anwm+1
                     IF(indwm(i1,i2,i3+1)) anwm=anwm+1
                     IF(anwm.ge.n0) THEN
C reassign to white matter
                        anz=anz+1
                        indwm(i1,i2,i3) = .TRUE.
                        indre(i1,i2,i3) = .FALSE.
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine detctbrd(simg,n1,n2,n3,dw,bimg)
      implicit logical (a-z)
      integer n1,n2,n3,simg(n1,n2,n3),dw
      logical bimg(n1,n2,n3),border
      integer i1,i2,i3,is
      DO i1=1,n1
         DO i2=1,n2
           DO i3=1,n3
              bimg(i1,i2,i3)=.FALSE.
           END DO
         END DO
      END DO
      DO i1=dw+1,n1-dw
         DO i2=dw+1,n2-dw
            DO i3=dw+1,n3-dw
               border=.FALSE.
               is=simg(i1,i2,i3)
               if(simg(i1-1,i2,i3).ne.is) border=.TRUE.
               if(simg(i1+1,i2,i3).ne.is) border=.TRUE.
               if(simg(i1,i2-1,i3).ne.is) border=.TRUE.
               if(simg(i1,i2+1,i3).ne.is) border=.TRUE.
               if(simg(i1,i2,i3-1).ne.is) border=.TRUE.
               if(simg(i1,i2,i3+1).ne.is) border=.TRUE.
               bimg(i1,i2,i3)=border
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine detctbr0(simg,n1,n2,n3,dw,bimg)
      implicit logical (a-z)
      integer n1,n2,n3,simg(n1,n2,n3),dw
      logical bimg(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3,is
      DO i1=1,n1
         DO i2=1,n2
           DO i3=1,n3
              bimg(i1,i2,i3)=.FALSE.
           END DO
         END DO
      END DO
      DO i1=dw+1,n1-dw
         DO i2=dw+1,n2-dw
            DO i3=dw+1,n3-dw
               is=simg(i1,i2,i3)
               outerloop: DO j1=i1-dw,i1+dw
                  DO j2=i2-dw,i2+dw
                     DO j3=i3-dw,i3+dw
                        if(simg(j1,j2,j3).ne.is) THEN
                           bimg(i1,i2,i3)=.TRUE.
                           exit outerloop
                        END IF
                     END DO
                  END DO
               END DO outerloop
            END DO
         END DO
      END DO
      RETURN
      END