
c
c
c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  wave,s,amdq,apdq)
c     =====================================================
c
c     # Roe-solver for the Euler equations with a tracer variable
c     # and separate shear and entropy waves.
c
c     # solve Riemann problems along one slice of data.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # This data is along a slice in the x-direction if ixy=1 
c     #                            or the y-direction if ixy=2.
c     # On output, wave contains the waves, s the speeds, 
c     # and amdq, apdq the decomposition of the flux difference
c     #   f(qr(i-1)) - f(ql(i))  
c     # into leftgoing and rightgoing parts respectively.
c     # With the Roe solver we have   
c     #    amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
c     # where A is the Roe matrix.  An entropy fix can also be incorporated
c     # into the flux differences.
c
c     # Note that the i'th Riemann problem has left state qr(:,i-1)
c     #                                    and right state ql(:,i)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit double precision (a-h,o-z)
c
      dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
      dimension    s(mwaves, 1-mbc:maxm+mbc)
      dimension   ql(meqn, 1-mbc:maxm+mbc)
      dimension   qr(meqn, 1-mbc:maxm+mbc)
      dimension  apdq(meqn, 1-mbc:maxm+mbc)
      dimension  amdq(meqn, 1-mbc:maxm+mbc)
c
c     local arrays -- common block comroe is passed to rpt2eu
c     ------------
      parameter (maxm2 = 1800)  
      dimension delta(4)
      logical efix
      common /cparam/  gamma,gamma1
c     # assumes at most maxm2 * maxm2 grid with mbc<=3
c
      data efix /.true./    !# use entropy fix for transonic rarefactions
c
c
c     # set mu to point to  the component of the system that corresponds
c     # to momentum in the direction of this slice, mv to the orthogonal
c     # momentum:
c
      if (ixy.eq.1) then
          mu = 2
          mv = 3
        else
          mu = 3
          mv = 2
        endif
c
c     # note that notation for u and v reflects assumption that the 
c     # Riemann problems are in the x-direction with u in the normal
c     # direciton and v in the orthogonal direcion, but with the above
c     # definitions of mu and mv the routine also works with ixy=2
c     # and returns, for example, f0 as the Godunov flux g0 for the
c     # Riemann problems u_t + g(u)_y = 0 in the y-direction.
c
c
c     # compute the Roe-averaged variables needed in the Roe solver.
c     # These are stored in the common block comroe since they are
c     # later used in routine rpt2eu to do the transverse wave splitting.
c
      do 20 i = 2-mbc, mx+mbc
         rhsqrtl = dsqrt(qr(1,i-1))
         rhsqrtr = dsqrt(ql(1,i))
         pl = gamma1*(qr(4,i-1) - 0.5d0*(qr(2,i-1)**2 +
     &        qr(3,i-1)**2)/qr(1,i-1))
         pr = gamma1*(ql(4,i) - 0.5d0*(ql(2,i)**2 +
     &        ql(3,i)**2)/ql(1,i))
         rhsq2 = rhsqrtl + rhsqrtr
         u = (qr(mu,i-1)/rhsqrtl + ql(mu,i)/rhsqrtr) / rhsq2
         v = (qr(mv,i-1)/rhsqrtl + ql(mv,i)/rhsqrtr) / rhsq2
         enth = (((qr(4,i-1)+pl)/rhsqrtl
     &             + (ql(4,i)+pr)/rhsqrtr)) / rhsq2
         u2v2 = u**2 + v**2
         a2 = gamma1*(enth - .5d0*u2v2)
         a = dsqrt(a2)
         g1a2 = gamma1 / a2
         euv = enth - u2v2 
c
c
c     # now split the jump in q at each interface into waves
c
c     # find a1 thru a4, the coefficients of the 4 eigenvectors:
         delta(1) = ql(1,i) - qr(1,i-1)
         delta(2) = ql(mu,i) - qr(mu,i-1)
         delta(3) = ql(mv,i) - qr(mv,i-1)
         delta(4) = ql(4,i) - qr(4,i-1)
         a3 = g1a2 * (euv*delta(1) 
     &      + u*delta(2) + v*delta(3) - delta(4))
         a2 = delta(3) - v*delta(1)
         a4 = (delta(2) + (a-u)*delta(1) - a*a3) / (2.d0*a)
         a1 = delta(1) - a3 - a4
c
c        # Compute the waves.
c
c        # acoustic:
         wave(1,1,i) = a1
         wave(mu,1,i) = a1*(u-a)
         wave(mv,1,i) = a1*v
         wave(4,1,i) = a1*(enth - u*a)
         wave(5,1,i) = 0.d0
         s(1,i) = u-a
c
c        # shear:
         wave(1,2,i) = 0.d0
         wave(mu,2,i) = 0.d0
         wave(mv,2,i) = a2
         wave(4,2,i) = a2*v
         wave(5,2,i) = 0.d0
         s(2,i) = u
c
c        # entropy:
         wave(1,3,i) = a3
         wave(mu,3,i) = a3*u
         wave(mv,3,i) = a3*v
         wave(4,3,i) = a3*0.5d0*u2v2
         wave(5,3,i) = 0.d0
         s(3,i) = u
c
c        # acoustic:
         wave(1,4,i) = a4
         wave(mu,4,i) = a4*(u+a)
         wave(mv,4,i) = a4*v
         wave(4,4,i) = a4*(enth+u*a)
         wave(5,4,i) = 0.d0
         s(4,i) = u+a
c
c        # Another wave added for tracer concentration:
c
c        # tracer:
         wave(1,5,i) = 0.d0
         wave(mu,5,i) = 0.d0
         wave(mv,5,i) = 0.d0
         wave(4,5,i) = 0.d0
         wave(5,5,i) = ql(5,i) - qr(5,i-1)
         s(5,i) = u
c
   20    continue
c
c
c    # compute flux differences amdq and apdq.
c    ---------------------------------------
c
      if (efix) go to 110
c
c     # no entropy fix
c     ----------------
c
c     # amdq = SUM s*wave   over left-going waves
c     # apdq = SUM s*wave   over right-going waves
c
      do 100 m=1,meqn
         do 100 i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do 90 mw=1,mwaves
               if (s(mw,i) .lt. 0.d0) then
                   amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                 else
                   apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                 endif
   90          continue
  100       continue
      go to 900     
c
c-----------------------------------------------------
c
  110 continue
c
c     # With entropy fix
c     ------------------
c
c    # compute flux differences amdq and apdq.
c    # First compute amdq as sum of s*wave for left going waves.
c    # Incorporate entropy fix by adding a modified fraction of wave
c    # if s should change sign.
c
      do 200 i = 2-mbc, mx+mbc
c
c        # check 1-wave:
c        ---------------
c
         rhoim1 = qr(1,i-1)
         pim1 = gamma1*(qr(4,i-1) - 0.5d0*(qr(mu,i-1)**2 
     &           + qr(mv,i-1)**2) / rhoim1)
         cim1 = dsqrt(gamma*pim1/rhoim1)
         s0 = qr(mu,i-1)/rhoim1 - cim1     !# u-c in left state (cell i-1)

c        # check for fully supersonic case:
         if (s0.ge.0.d0 .and. s(1,i).gt.0.d0)  then
c            # everything is right-going
             do 60 m=1,meqn
                amdq(m,i) = 0.d0
   60           continue
             go to 200 
             endif
c
         rho1 = qr(1,i-1) + wave(1,1,i)
         rhou1 = qr(mu,i-1) + wave(mu,1,i)
         rhov1 = qr(mv,i-1) + wave(mv,1,i)
         en1 = qr(4,i-1) + wave(4,1,i)
         p1 = gamma1*(en1 - 0.5d0*(rhou1**2 + rhov1**2)/rho1)
         c1 = dsqrt(gamma*p1/rho1)
         s1 = rhou1/rho1 - c1  !# u-c to right of 1-wave
         if (s0.lt.0.d0 .and. s1.gt.0.d0) then
c            # transonic rarefaction in the 1-wave
             sfract = s0 * (s1-s(1,i)) / (s1-s0)
           else if (s(1,i) .lt. 0.d0) then
c            # 1-wave is leftgoing
             sfract = s(1,i)
           else
c            # 1-wave is rightgoing
             sfract = 0.d0   !# this shouldn't happen since s0 < 0
           endif
         do 120 m=1,meqn
            amdq(m,i) = sfract*wave(m,1,i)
  120       continue
c
c        # check contact discontinuity:
c        ------------------------------
c
         if (s(2,i) .ge. 0.d0) go to 200  !# 2- 3- and 5-waves are rightgoing
         do 140 m=1,meqn
            amdq(m,i) = amdq(m,i) + s(2,i)*wave(m,2,i)
            amdq(m,i) = amdq(m,i) + s(3,i)*wave(m,3,i)
            amdq(m,i) = amdq(m,i) + s(5,i)*wave(m,5,i)
  140       continue
c
c        # check 4-wave:
c        ---------------
c
         rhoi = ql(1,i)
         pi = gamma1*(ql(4,i) - 0.5d0*(ql(mu,i)**2 
     &           + ql(mv,i)**2) / rhoi)
         ci = dsqrt(gamma*pi/rhoi)
         s3 = ql(mu,i)/rhoi + ci     !# u+c in right state  (cell i)
c
         rho2 = ql(1,i) - wave(1,4,i)
         rhou2 = ql(mu,i) - wave(mu,4,i)
         rhov2 = ql(mv,i) - wave(mv,4,i)
         en2 = ql(4,i) - wave(4,4,i)
         p2 = gamma1*(en2 - 0.5d0*(rhou2**2 + rhov2**2)/rho2)
         c2 = dsqrt(gamma*p2/rho2)
         s2 = rhou2/rho2 + c2   !# u+c to left of 4-wave
         if (s2 .lt. 0.d0 .and. s3.gt.0.d0) then
c            # transonic rarefaction in the 4-wave
             sfract = s2 * (s3-s(4,i)) / (s3-s2)
           else if (s(4,i) .lt. 0.d0) then
c            # 4-wave is leftgoing
             sfract = s(4,i)
           else 
c            # 4-wave is rightgoing
             go to 200
           endif
c
         do 160 m=1,meqn
            amdq(m,i) = amdq(m,i) + sfract*wave(m,4,i)
  160       continue
  200    continue
c
c     # compute the rightgoing flux differences:
c     # df = SUM s*wave   is the total flux difference and apdq = df - amdq
c
      do 220 m=1,meqn
         do 220 i = 2-mbc, mx+mbc
            df = 0.d0
            do 210 mw=1,mwaves
               df = df + s(mw,i)*wave(m,mw,i)
  210          continue
            apdq(m,i) = df - amdq(m,i)
  220       continue
c
  900 continue
      return
      end
