c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  imp,asdq,bmasdq,bpasdq)
c     =====================================================
c
c     # Riemann solver in the transverse direction.
c     # This is a dummy routine that returns zeros and is only intended
c     # to illustrate the format of this routine.  See various example
c     # directories for better examples.

c     # This dummy routine can be used if transverse solves are not being
c     # used, i.e. if method(3)=0.
c
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c
      implicit double precision (a-h,o-z)
      dimension     ql(meqn,1-mbc:maxm+mbc)
      dimension     qr(meqn,1-mbc:maxm+mbc)
      dimension   asdq(meqn,1-mbc:maxm+mbc)
      dimension bmasdq(meqn,1-mbc:maxm+mbc)
      dimension bpasdq(meqn,1-mbc:maxm+mbc)

      do 10 i = 2-mbc, mx+mbc
         !down going part of the normal fluctuation
         bmasdq(1,i)=0.d0
         bmasdq(2,i)=0.d0
         bmasdq(3,i)=0.d0
         !up going part of the normal fluctuation
         bpasdq(1,i)=0.d0
         bpasdq(2,i)=0.d0
         bpasdq(3,i)=0.d0
 10      continue

      return
      end
