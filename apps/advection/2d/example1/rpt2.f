c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  ilr,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Riemann solver in the transverse direction for the scalar equation
c
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c
      dimension     ql(meqn,1-mbc:maxm+mbc)
      dimension     qr(meqn,1-mbc:maxm+mbc)
      dimension   asdq(meqn,1-mbc:maxm+mbc)
      dimension bmasdq(meqn,1-mbc:maxm+mbc)
      dimension bpasdq(meqn,1-mbc:maxm+mbc)
      parameter (maxm2 = 502)
      common /comrp/ ubar,vbar
c
c     # transverse wave speeds have been computed in rpn2
c     # maux=0 and aux arrays are unused in this example.
c
      if (ixy.eq.1) then
             stran = vbar
           else
             stran = ubar
           endif

      stranm = dmin1(stran, 0.d0)
      stranp = dmax1(stran, 0.d0)
      
      do 10 i = 2-mbc, mx+mbc
          bmasdq(1,i) = stranm * asdq(1,i)
          bpasdq(1,i) = stranp * asdq(1,i)
   10     continue
c
      return
      end
