
c
c
c
c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,
     &                  ylower,dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
       implicit double precision (a-h,o-z)
       dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
       dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)
       common /dmetric/ dmetric(-1:300+2,-1:150+2,5)

c
       pi2 = 8.d0*datan(1.d0)  !# = 2 * pi
       pi = 0.5*pi2

         do 20 i=1-mbc,mx+mbc
	   xc = xlower + (i-0.5d0)*dx
           do 20 j=1-mbc,my+mbc
	      yc = ylower + (j-0.5d0)*dy
	      call mapc2m(xc,yc,xp,yp,zp)
	      if (xp.gt.0.d0 .and. yp.gt.0.d0) then
                 q(i,j,1) = 1.d0
		else
                 q(i,j,1) = 0.d0
		endif
  20   continue
 

       return
       end
