c =========================================================
      subroutine rp1(maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &               wave,s,amdq,apdq)
c =========================================================
c
c     # Solve Riemann problems for the 1D shallow water equations
c     #   (h)_t + (u h)_x = 0 
c     #   (uh)_t + ( uuh + .5*gh^2 )_x = 0 
c     # using Roe's approximate Riemann solver with entropy fix for
c     # transonic rarefractions.  
c     
c     # This function solve the Riemann problem at all interfaces in one
c     # call
c
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
c     # Here meqn=mwaves=2 should be passed from the calling routine
c
c     # This Riemann solver differs from the original clawpack Riemann solver 
c     # for the interleaved indices



      implicit double precision (a-h,o-z)
      dimension   ql(meqn,         1-mbc:maxmx+mbc)
      dimension   qr(meqn,         1-mbc:maxmx+mbc)
      dimension    s(mwaves,       1-mbc:maxmx+mbc)
      dimension wave(meqn, mwaves, 1-mbc:maxmx+mbc)
      dimension amdq(meqn,         1-mbc:maxmx+mbc)
      dimension apdq(meqn,         1-mbc:maxmx+mbc)

c     # Gravity constant set in the shallow1D.py file      
      common /cparam/ grav

c
c     # Local storage
c     ---------------
      dimension delta(2)
      logical efix
c     common /comrp/ grav
      data efix /.true./    !# use entropy fix for transonic rarefactions




c     # Main loop of the Riemann solver. 
      do 30 i=2-mbc,mx+mbc
c
c         
c     # compute  Roe-averaged quantities: 
         ubar = (qr(2,i-1)/dsqrt(qr(1,i-1)) + ql(2,i)/dsqrt(ql(1,i)))/
     .       ( dsqrt(qr(1,i-1)) + dsqrt(ql(1,i)) )
         cbar=dsqrt(0.5d0*grav*(qr(1,i-1) + ql(1,i)))
         
c     # delta(1)=h(i)-h(i-1) and  delta(2)=hu(i)-hu(i-1)
      delta(1) = ql(1,i) - qr(1,i-1)
      delta(2) = ql(2,i) - qr(2,i-1)

c     # Compute coeffs in the evector expansion of delta(1),delta(2)
      a1 = 0.5d0*(-delta(2) + (ubar + cbar) * delta(1))/cbar
      a2 = 0.5d0*( delta(2) - (ubar - cbar) * delta(1))/cbar

c     # Finally, compute the waves.
      wave(1,1,i) = a1
      wave(2,1,i) = a1*(ubar - cbar)
      s(1,i) = ubar - cbar
         
      wave(1,2,i) = a2
      wave(2,2,i) = a2*(ubar + cbar)
      s(2,i) = ubar + cbar
         
   30 continue

c     # Compute Godunov flux f0:
c     --------------------------

      if (efix) go to 110

c     # No entropy fix
c     ----------------------------------------------
c     # amdq = SUM s*wave   over left-going waves
c     # apdq = SUM s*wave   over right-going waves

      do 100 m=1,2
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
c    -----------------------------------------------
c

  110 continue

c     # With entropy fix
c     ------------------
c
c    # compute flux differences amdq and apdq.
c    # First compute amdq as sum of s*wave for left going waves.
c    # Incorporate entropy fix by adding a modified fraction of wave
c    # if s should change sign.
c
      do 200 i=2-mbc,mx+mbc
      
c ------------------------------------------------------
c        # check 1-wave:
c        ---------------
c
c        # u-c in left state (cell i-1)
         s0 = qr(2,i-1)/qr(1,i-1) - dsqrt(grav*qr(1,i-1))
         
c        # check for fully supersonic case:
         if (s0.ge.0.d0 .and. s(1,i).gt.0.d0)  then
c            # everything is right-going
             do 60 m=1,2
                amdq(m,i) = 0.d0
   60        continue
             go to 200 
         endif
c
c        # u-c to right of 1-wave
         hr1  = qr(1,i-1) + wave(1,1,i) 
         uhr1 = qr(2,i-1) + wave(2,1,i)
         s1 =  uhr1/hr1 - dsqrt(grav*hr1)
                 
         if (s0.lt.0.d0 .and. s1.gt.0.d0) then
c            # transonic rarefaction in the 1-wave
             sfract = s0 * (s1-s(1,i)) / (s1-s0)
         else if (s(1,i) .lt. 0.d0) then
c	     # 1-wave is leftgoing
             sfract = s(1,i)
         else
c	     # 1-wave is rightgoing
             sfract = 0.d0   !# this shouldn't happen since s0 < 0
         endif

         do 120 m=1,2
             amdq(m,i) = sfract*wave(m,1,i)
  120    continue
  
c -------------------------------------------------------
c        # check 2-wave:
c        ---------------
c        # u+c in right state  (cell i)
         s3 = ql(2,i)/ql(1,i) + dsqrt(grav*ql(1,i))
              
c        # u+c to left of 2-wave
         hl2  = ql(1,i) - wave(1,2,i) 
         uhl2 = ql(2,i) - wave(2,2,i)
         s2 = uhl2/hl2 + dsqrt(grav*hl2)
                  
         if (s2 .lt. 0.d0 .and. s3.gt.0.d0) then
c            # transonic rarefaction in the 2-wave
             sfract = s2 * (s3-s(2,i)) / (s3-s2)
         else if (s(2,i) .lt. 0.d0) then
c            # 2-wave is leftgoing
             sfract = s(2,i)
         else 
c            # 2-wave is rightgoing
             go to 200
        endif
c
        do 160 m=1,2
            amdq(m,i) = amdq(m,i) + sfract*wave(m,2,i)
  160       continue
  200    continue
c
c
c     # compute the rightgoing flux differences:
c     # df = SUM s*wave   is the total flux difference and apdq = df - amdq
c
      do 220 m=1,2
         do 220 i = 2-mbc, mx+mbc
            df = 0.d0
            do 210 mw=1,mwaves
               df = df + s(mw,i)*wave(m,mw,i)
  210          continue
            apdq(m,i) = df - amdq(m,i)
  220       continue

  900 continue
      return
      end



