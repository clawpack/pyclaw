c     =====================================================
      subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &     dx,dy,q,maux,aux)
c     =====================================================
c
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      common /cqinit/ A1,beta1,x1,y1, A2,beta2,x2,y2

      do i=1,mx
        xc = xlower + (i-0.5d0)*dx
        do j=1,my
            yc = ylower + (j-0.5d0)*dy
            call mapc2p(xc,yc,xp,yp)
            q(i,j,1) = A1*exp(-beta1*((xp - x1)**2 + (yp - y1)**2)) 
     &                  +A2*exp(-beta2*((xp - x2)**2 + (yp - y2)**2))
        enddo
      enddo

      return
      end
