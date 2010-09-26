c
c
c =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c     # Smooth entropy wave hitting a shock
c
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)
      common /param/ gamma, gamma1
c
c
      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx
         q(i,1) = 1.d0 + 0.5d0*dexp(-80.d0*xcell**2)
         q(i,2) = 0.d0
         q(i,3) = q(i,1) / gamma1
  150    continue
c
      return
      end
