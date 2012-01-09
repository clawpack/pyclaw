c
c
c     =====================================================
      subroutine limiter(maxm,meqn,mwaves,mbc,mx,wave,s,mthlim)
c     =====================================================
c
c     # Apply a limiter to the waves.
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
c     # value of the corresponding element of the array mthlim, as used in
c     # the function philim.
c     # Note that a different limiter may be used in each wave family.
c
c     # dotl and dotr denote the inner product of wave with the wave to
c     # the left or right.  The norm of the projections onto the wave are then
c     # given by dotl/wnorm2 and dotr/wnorm2, where wnorm2 is the 2-norm
c     # of wave.
c
      implicit double precision (a-h,o-z)
      dimension mthlim(mwaves)
      dimension wave(meqn, mwaves,1-mbc:maxm+mbc)
      dimension    s(mwaves,1-mbc:maxm+mbc)
c
c
      do 50 mw=1,mwaves
         if (mthlim(mw) .eq. 0) go to 50
         dotr = 0.d0
         do 40 i = 0, mx+1
            wnorm2 = 0.d0
            dotl = dotr
            dotr = 0.d0
            do 20 m=1,meqn
               wnorm2 = wnorm2 + wave(m,mw,i)**2
               dotr = dotr + wave(m,mw,i)*wave(m,mw,i+1)
   20          continue
            if (i.eq.0) go to 40
            if (wnorm2.eq.0.d0) go to 40
c
            if (s(mw,i) .gt. 0.d0) then
                wlimitr = philim(wnorm2, dotl, mthlim(mw))
              else
                wlimitr = philim(wnorm2, dotr, mthlim(mw))
              endif
c
            do 30 m=1,meqn
               wave(m,mw,i) = wlimitr * wave(m,mw,i)
   30          continue
   40       continue
   50    continue
c
      return
      end
