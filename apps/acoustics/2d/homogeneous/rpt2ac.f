c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  imp,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Riemann solver in the transverse direction for the acoustics equations.
c
c     # Split asdq into down-going flux bmasdq and up-going flux bpasdq.
c
      dimension     ql(meqn, 1-mbc:maxm+mbc)
      dimension     qr(meqn, 1-mbc:maxm+mbc)
      dimension   asdq(meqn, 1-mbc:maxm+mbc)
      dimension bmasdq(meqn, 1-mbc:maxm+mbc)
      dimension bpasdq(meqn, 1-mbc:maxm+mbc)
c
c     # density, bulk modulus, and sound speed, and impedence of medium:
c     # (should be set in setprob.f)
      common /cparam/ rho,bulk,cc,zz

c
      if (ixy.eq.1) then
          mu = 2
          mv = 3
        else
          mu = 3
          mv = 2
        endif
c
      do 20 i = 2-mbc, mx+mbc
         a1 = (-asdq(1,i) + zz*asdq(mv,i)) / (2.d0*zz)
         a2 = (asdq(1,i) + zz*asdq(mv,i)) / (2.d0*zz)
c
c        # The down-going flux difference bmasdq is the product  -c * wave
c
         bmasdq(1,i) = cc * a1*zz
         bmasdq(mu,i) = 0.d0
         bmasdq(mv,i) = -cc * a1
c
c        # The up-going flux difference bpasdq is the product  c * wave
c
         bpasdq(1,i) = cc * a2*zz
         bpasdq(mu,i) = 0.d0
         bpasdq(mv,i) = cc * a2
c
   20    continue
c
      return
      end
