c
c
c =========================================================
      subroutine rp1(maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &           wave,s,amdq,apdq)
c =========================================================
c
c     # solve Riemann problems for the 1D advection equation q_t + u*q_x = 0.
c     # For constant advection velocity u, passed in common block.
c
c     # The advection speed u is passed in the common block comrp
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c     # On output, wave contains the waves,
c     #            s the speeds,
c     #            amdq the  left-going flux difference  A^- \Delta q
c     #            apdq the right-going flux difference  A^+ \Delta q
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routine step1, rp is called with ql = qr = q.
c
c
      implicit double precision (a-h,o-z)
      dimension   ql( meqn,1-mbc:maxmx+mbc)
      dimension   qr( meqn,1-mbc:maxmx+mbc)
      dimension   amdq( meqn,1-mbc:maxmx+mbc)
      dimension   apdq( meqn,1-mbc:maxmx+mbc)
      dimension    s(mwaves,1-mbc:maxmx+mbc)
      dimension wave(meqn, mwaves,1-mbc:maxmx+mbc)
      common /comrp/ u
c
c
      do 30 i=2-mbc,mx+mbc
c
c        # Compute the wave and speed
c
         wave(1,1,i) = ql(1,i) - qr(1,i-1)
         s(1,i) = u
         amdq(1,i) = dmin1(u, 0.d0) * wave(1,1,i)
         apdq(1,i) = dmax1(u, 0.d0) * wave(1,1,i)
   30    continue
c
      return
      end
