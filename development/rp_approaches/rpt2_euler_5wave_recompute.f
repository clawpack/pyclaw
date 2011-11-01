c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  ilr,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Riemann solver in the transverse direction for the Euler equations
c     #  with a tracer variable.
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c     # Uses Roe averages and other quantities which were 
c     # computed in rpn2eu and stored in the common block comroe.
c
      dimension     ql(meqn, 1-mbc:maxm+mbc)
      dimension     qr(meqn, 1-mbc:maxm+mbc)
      dimension   asdq(meqn, 1-mbc:maxm+mbc)
      dimension bmasdq(meqn, 1-mbc:maxm+mbc)
      dimension bpasdq(meqn, 1-mbc:maxm+mbc)
c
      common /cparam/  gamma,gamma1
      dimension waveb(5,4),sb(4)
      parameter (maxm2 = 1800)  
c     # assumes at most 800x600 grid with mbc<=3
      dimension u2v2(-2:maxm2+3),
     &       u(-2:maxm2+3),v(-2:maxm2+3),
     &       enth(-2:maxm2+3),a(-2:maxm2+3),
     &       g1a2(-2:maxm2+3),euv(-2:maxm2+3) 
c
      if (mbc.gt.3 .or. maxm2 .lt. maxm) then
         write(6,*) 'need to increase maxm2 or 3 in rpt'
         stop
         endif
c
      if (ixy.eq.1) then
          mu = 2
          mv = 3
        else
          mu = 3
          mv = 2
        endif


c     # compute the Roe-averaged variables needed in the Roe solver.
c
      do 10 i = 2-mbc, mx+mbc
         rhsqrtl = dsqrt(qr(1,i-1))
         rhsqrtr = dsqrt(ql(1,i))
         pl = gamma1*(qr(4,i-1) - 0.5d0*(qr(2,i-1)**2 +
     &        qr(3,i-1)**2)/qr(1,i-1))
         pr = gamma1*(ql(4,i) - 0.5d0*(ql(2,i)**2 +
     &        ql(3,i)**2)/ql(1,i))
         rhsq2 = rhsqrtl + rhsqrtr
         u(i) = (qr(mu,i-1)/rhsqrtl + ql(mu,i)/rhsqrtr) / rhsq2
         v(i) = (qr(mv,i-1)/rhsqrtl + ql(mv,i)/rhsqrtr) / rhsq2
         enth(i) = (((qr(4,i-1)+pl)/rhsqrtl
     &             + (ql(4,i)+pr)/rhsqrtr)) / rhsq2
         u2v2(i) = u(i)**2 + v(i)**2
         a2 = gamma1*(enth(i) - .5d0*u2v2(i))
         a(i) = dsqrt(a2)
         g1a2(i) = gamma1 / a2
         euv(i) = enth(i) - u2v2(i) 
   10    continue
c
c
         do 20 i = 2-mbc, mx+mbc
            a3 = g1a2(i) * (euv(i)*asdq(1,i) 
     &             + u(i)*asdq(mu,i) + v(i)*asdq(mv,i) - asdq(4,i))
            a2 = asdq(mu,i) - u(i)*asdq(1,i)
            a4 = (asdq(mv,i) + (a(i)-v(i))*asdq(1,i) - a(i)*a3)
     &              / (2.d0*a(i))
            a1 = asdq(1,i) - a3 - a4
c
            waveb(1,1) = a1
            waveb(mu,1) = a1*u(i)
            waveb(mv,1) = a1*(v(i)-a(i))
            waveb(4,1) = a1*(enth(i) - v(i)*a(i))
            waveb(5,1) = 0.d0
            sb(1) = v(i) - a(i)
c
            waveb(1,2) = a3
            waveb(mu,2) = a3*u(i) + a2
            waveb(mv,2) = a3*v(i)
            waveb(4,2) = a3*0.5d0*u2v2(i) + a2*u(i)
            waveb(5,2) = 0.d0
            sb(2) = v(i)
c
            waveb(1,3) = a4
            waveb(mu,3) = a4*u(i)
            waveb(mv,3) = a4*(v(i)+a(i))
            waveb(4,3) = a4*(enth(i)+v(i)*a(i))
            waveb(5,3) = 0.d0
            sb(3) = v(i) + a(i)
c
            waveb(1,4) = 0.d0
            waveb(mu,4) = 0.d0
            waveb(mv,4) = 0.d0
            waveb(4,4) = 0.d0
            waveb(5,4) = asdq(5,i)
            sb(4) = v(i)
c
c           # compute the flux differences bmasdq and bpasdq
c
            do 30 m=1,meqn
               bmasdq(m,i) = 0.d0
               bpasdq(m,i) = 0.d0
               do 30 mw=1,4
                  bmasdq(m,i) = bmasdq(m,i) 
     &                         + dmin1(sb(mw), 0.d0) * waveb(m,mw)
                  bpasdq(m,i) = bpasdq(m,i)
     &                         + dmax1(sb(mw), 0.d0) * waveb(m,mw)
   30             continue
c                 
   20       continue
c
      return
      end
