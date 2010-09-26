c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  imp,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Riemann solver in the transverse direction for 2D Burgers' equation
c
c     # Split asdq into eigenvectors of Roe matrix B.
c     # For the scalar equation, this simply amounts to computing the
c     # transverse wave speed from the opposite Riemann problem.
c
      dimension    ql(1-mbc:maxm+mbc, meqn)
      dimension    qr(1-mbc:maxm+mbc, meqn)
      dimension   asdq(1-mbc:maxm+mbc, meqn)
      dimension bmasdq(1-mbc:maxm+mbc, meqn)
      dimension bpasdq(1-mbc:maxm+mbc, meqn)
c
c     # x- and y- Riemann problems are identical, so it doesn't matter if
c     # ixy=1 or 2.
c
          do 10 i = 2-mbc, mx+mbc
             sb = 0.5d0*(qr(i-1,1) + ql(i,1))
             bmasdq(i,1) = dmin1(sb, 0.d0) * asdq(i,1)
             bpasdq(i,1) = dmax1(sb, 0.d0) * asdq(i,1)
   10        continue
c
      return
      end
