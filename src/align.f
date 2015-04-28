      subroutine align3d(Ms,Rsinv,Rt,off,nxs,nys,nzs,nxt,nyt,nzt,Mt)
C
C  Realign without interpolation
C
      implicit logical (a-z)
      integer nxs,nxt,nys,nyt,nzs,nzt
      real*8 Mt(nxt,nyt,nzt),Ms(nxs,nys,nzs),Rsinv(3,3),Rt(3,3),off(3)
      integer tx,ty,tz,sx,sy,sz
      real*8 zx,zy,zz
      DO tx=0,nxt-1 
         DO ty=0,nyt-1
            DO tz=0,nzt-1
               Mt(tx+1,ty+1,tz+1)=0.d0
               zx=Rt(1,1)*tx+Rt(1,2)*ty+Rt(1,3)*tz+off(1)
               zy=Rt(2,1)*tx+Rt(2,2)*ty+Rt(2,3)*tz+off(2)
               zz=Rt(3,1)*tx+Rt(3,2)*ty+Rt(3,3)*tz+off(3)
               sx=Rsinv(1,1)*zx+Rsinv(1,2)*zy+Rsinv(1,3)*zz+1.d0
               sy=Rsinv(2,1)*zx+Rsinv(2,2)*zy+Rsinv(2,3)*zz+1.d0
               sz=Rsinv(3,1)*zx+Rsinv(3,2)*zy+Rsinv(3,3)*zz+1.d0
               IF(sx.lt.1.or.sx.gt.nxs) CYCLE
               IF(sy.lt.1.or.sy.gt.nys) CYCLE
               IF(sz.lt.1.or.sz.gt.nzs) CYCLE
               Mt(tx+1,ty+1,tz+1)=Ms(sx,sy,sz)
            END DO
         END DO
      END DO
      RETURN
      END