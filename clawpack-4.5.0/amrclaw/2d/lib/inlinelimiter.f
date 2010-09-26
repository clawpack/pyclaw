c
c
c     =====================================================
      subroutine limiter(maxm,meqn,mwaves,mbc,mx,wave,s,mthlim)
c     =====================================================
c
c     # Apply a limiter to the waves.  
c
c     # Version of December, 2002.
c     # Modified from the original CLAWPACK routine to eliminate calls 
c     # to philim.  Since philim was called for every wave at each cell
c     # interface, this was adding substantial overhead in some cases.
c
c     # The limiter is computed by comparing the 2-norm of each wave with
c     # the projection of the wave from the interface to the left or
c     # right onto the current wave.  For a linear system this would
c     # correspond to comparing the norms of the two waves.  For a 
c     # nonlinear problem the eigenvectors are not colinear and so the 
c     # projection is needed to provide more limiting in the case where the
c     # neighboring wave has large norm but points in a different direction
c     # in phase space.
c
c     # The specific limiter used in each family is determined by the
c     # value of the corresponding element of the array mthlim.
c     # Note that a different limiter may be used in each wave family.
c
c     # dotl and dotr denote the inner product of wave with the wave to
c     # the left or right.  The norm of the projections onto the wave are then
c     # given by dotl/wnorm2 and dotr/wnorm2, where wnorm2 is the 2-norm
c     # of wave.
c
      implicit double precision (a-h,o-z)
      dimension mthlim(mwaves)
      dimension wave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
c
c
      do 200 mw=1,mwaves
         if (mthlim(mw) .eq. 0) go to 200
         dotr = 0.d0
         do 190 i = 0, mx+1
            wnorm2 = 0.d0
            dotl = dotr
            dotr = 0.d0
            do 5 m=1,meqn
               wnorm2 = wnorm2 + wave(i,m,mw)**2
               dotr = dotr + wave(i,m,mw)*wave(i+1,m,mw)
    5          continue
            if (i.eq.0) go to 190
            if (wnorm2.eq.0.d0) go to 190
c
            if (s(i,mw) .gt. 0.d0) then
                r = dotl / wnorm2
              else
                r = dotr / wnorm2
              endif
c
            go to (10,20,30,40,50) mthlim(mw)
c
   10       continue
c           --------
c           # minmod
c           --------
            wlimitr = dmax1(0.d0, dmin1(1.d0, r))
            go to 170
c
   20       continue
c           ----------
c           # superbee
c           ----------
            wlimitr = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
            go to 170
c
   30       continue
c           ----------
c           # van Leer
c           ----------
            wlimitr = (r + dabs(r)) / (1.d0 + dabs(r))
            go to 170
c
   40       continue
c           ------------------------------
c           # monotinized centered
c           ------------------------------
            c = (1.d0 + r)/2.d0
            wlimitr = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
            go to 170
c
   50       continue
c           ------------------------------
c           # Beam-Warming
c           ------------------------------
            wlimitr = r
            go to 170
c
  170       continue
c
c           # apply limiter to waves:
c
            do 180 m=1,meqn
               wave(i,m,mw) = wlimitr * wave(i,m,mw)
  180          continue

  190       continue
  200    continue
c
      return
      end
