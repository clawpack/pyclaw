c
c
c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                wave,s,amdq,apdq)
c     =====================================================
c
c     # Roe-solver for the 2D shallow water equations
c     # solve Riemann problems along one slice of data.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # This data is along a slice in the x-direction if ixy=1
c     #                            or the y-direction if ixy=2.
c
c     # This function solve the Riemann problem at all interfaces in one
c     # call
c
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
c     # This Riemann solver differs from the original clawpack Riemann solver 
c     # for the interleaved indices

      implicit double precision (a-h,o-z)

      dimension   ql(meqn,           1-mbc:maxm+mbc)
      dimension   qr(meqn,           1-mbc:maxm+mbc)
      dimension    s(mwaves,         1-mbc:maxm+mbc)
      dimension wave(meqn,   mwaves, 1-mbc:maxm+mbc)
      dimension  apdq(meqn,          1-mbc:maxm+mbc)
      dimension  amdq(meqn,          1-mbc:maxm+mbc)

c     # Gravity constant set in the shallow1D.py file 
      common /cparam/ grav

c     # Roe averages quantities of each interface
c     # These arrays are used afterwards in the transverse Riemann
c     # solver, i.e. rpt2sw.f
      common /comroe/ u(-2:603),v(-2:603),a(-2:603),h(-2:603)


c
c     local arrays
c     ------------
      dimension delta(3)
      logical efix

      data efix /.true./    !# Use entropy fix for transonic rarefactions

c     # Set mu to point to  the component of the system that corresponds
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
c     # Note that notation for u and v reflects assumption that the
c     # Riemann problems are in the x-direction with u in the normal
c     # direciton and v in the orthogonal direcion, but with the above
c     # definitions of mu and mv the routine also works with ixy=2
c     # and returns, for example, f0 as the Godunov flux g0 for the
c     # Riemann problems u_t + g(u)_y = 0 in the y-direction.
c
c
c     # Compute the Roe-averaged variables needed in the Roe solver.
c     # These are stored in the common block comroe since they are
c     # later used in routine rpt2sh to do the transverse wave splitting.
c
        do 10 i = 2-mbc, mx+mbc
         h(i) = (qr(1,i-1)+ql(1,i))*0.50d0
         hsqrtl = dsqrt(qr(1,i-1))
         hsqrtr = dsqrt(ql(1,i))
         hsq2 = hsqrtl + hsqrtr
         u(i) = (qr(mu,i-1)/hsqrtl + ql(mu,i)/hsqrtr) / hsq2
         v(i) = (qr(mv,i-1)/hsqrtl + ql(mv,i)/hsqrtr) / hsq2
         a(i) = dsqrt(grav*h(i))
   10    continue
c
c
c     # now split the jump in q at each interface into waves
c
c     # find a1 thru a3, the coefficients of the 3 eigenvectors:
      do 20 i = 2-mbc, mx+mbc
         delta(1) = ql(1,i) - qr(1,i-1)
         delta(2) = ql(mu,i) - qr(mu,i-1)
         delta(3) = ql(mv,i) - qr(mv,i-1)
         a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
         a2 = -v(i)*delta(1) + delta(3)
         a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))
c
c        # Compute the waves.
c
         wave(1,1,i) = a1
         wave(mu,1,i) = a1*(u(i)-a(i))
         wave(mv,1,i) = a1*v(i)
         s(1,i) = u(i)-a(i)
c
         wave(1,2,i) = 0.0d0
         wave(mu,2,i) = 0.0d0
         wave(mv,2,i) = a2
         s(2,i) = u(i)
c
         wave(1,3,i) = a3
         wave(mu,3,i) = a3*(u(i)+a(i))
         wave(mv,3,i) = a3*v(i)
         s(3,i) = u(i)+a(i)
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
      do 100 m=1,3
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
         do 200 i=2-mbc,mx+mbc
c           check 1-wave
            him1 = qr(1,i-1)
            s0 =  qr(mu,i-1)/him1 - dsqrt(grav*him1)
c           check for fully supersonic case :
            if (s0.gt.0.0d0.and.s(1,i).gt.0.0d0) then
               do 60 m=1,3
                  amdq(m,i)=0.0d0
   60          continue
               goto 200
            endif
c
            h1 = qr(1,i-1)+wave(1,1,i)
            hu1= qr(mu,i-1)+wave(mu,1,i)
            s1 = hu1/h1 - dsqrt(grav*h1) !speed just to right of 1-wave
            if (s0.lt.0.0d0.and.s1.gt.0.0d0) then
c              transonic rarefaction in 1-wave
               sfract = s0*((s1-s(1,i))/(s1-s0))
            else if (s(1,i).lt.0.0d0) then
c              1-wave is leftgoing
               sfract = s(1,i)
            else
c              1-wave is rightgoing
               sfract = 0.0d0
            endif
            do 120 m=1,3
               amdq(m,i) = sfract*wave(m,1,i)
  120       continue

c           check 2-wave
            if (s(2,i).gt.0.0d0) then
c              #2 and 3 waves are right-going
               go to 200
               endif

            do 140 m=1,3
               amdq(m,i) = amdq(m,i) + s(2,i)*wave(m,2,i)
  140       continue
c
c           check 3-wave
c
            hi = ql(1,i)
            s03 = ql(mu,i)/hi + dsqrt(grav*hi)
            h3=ql(1,i)-wave(1,3,i)
            hu3=ql(mu,i)-wave(mu,3,i)
            s3=hu3/h3 + dsqrt(grav*h3)
            if (s3.lt.0.0d0.and.s03.gt.0.0d0) then
c              transonic rarefaction in 3-wave
               sfract = s3*((s03-s(3,i))/(s03-s3))
            else if (s(3,i).lt.0.0d0) then
c              3-wave is leftgoing
               sfract = s(3,i)
            else
c              3-wave is rightgoing
               goto 200
            endif
            do 160 m=1,3
               amdq(m,i) = amdq(m,i) + sfract*wave(m,3,i)
  160       continue
  200       continue
c
c           compute rightgoing flux differences :
c
            do 220 m=1,3
               do 220 i = 2-mbc,mx+mbc
                  df = 0.0d0
                  do 210 mw=1,mwaves
                     df = df + s(mw,i)*wave(m,mw,i)
  210             continue
                  apdq(m,i)=df-amdq(m,i)
  220          continue
c
c
  900          continue
               return
               end


