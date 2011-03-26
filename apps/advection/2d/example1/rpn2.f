
c
c
c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann solver for the sample scalar equation
c     #  q_t + u*q_x + v*q_y = 0
c     
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # This data is along a slice in the x-direction if ixy=1 
c     #                            or the y-direction if ixy=2.
c     # On output, wave contains the waves, 
c     #            s the speeds, 
c     #
c     #            amdq = A^- Delta q, 
c     #            apdq = A^+ Delta q, 
c     #                   the decomposition of the flux difference
c     #                       f(qr(i-1)) - f(ql(i))
c     #                   into leftgoing and rightgoing parts respectively.
c     #               
c
c     # Note that the i'th Riemann problem has left state qr(:,i-1)
c     #                                    and right state ql(:,i)
c     # From the basic clawpack routines, this routine is called with ql = qr
c     # maux=0 and aux arrays are unused in this example.
c
c
      implicit double precision (a-h,o-z)
      common /comrp/ ubar,vbar
c
      dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
      dimension    s(mwaves, 1-mbc:maxm+mbc)
      dimension   ql(meqn, 1-mbc:maxm+mbc)
      dimension   qr(meqn, 1-mbc:maxm+mbc)
      dimension  apdq(meqn, 1-mbc:maxm+mbc)
      dimension  amdq(meqn, 1-mbc:maxm+mbc)
c
c
c$$$      do 30 i = 2-mbc, mx+mbc-1
      do 30 i = 2-mbc, mx+mbc
         wave(1,1,i) = ql(1,i) - qr(1,i-1)
         if (ixy.eq.1) then
             s(1,i) = ubar
           else
             s(1,i) = vbar
           endif
c
c        # flux differences:
         amdq(1,i) = dmin1(s(1,i), 0.d0) * wave(1,1,i)
         apdq(1,i) = dmax1(s(1,i), 0.d0) * wave(1,1,i)
c
   30    continue
c
      return
      end
