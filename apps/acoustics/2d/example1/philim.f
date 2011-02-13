c
c
c     =====================================================
      double precision function philim(a,b,meth)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Compute a limiter based on wave strengths a and b.
c     # meth determines what limiter is used.
c     # a is assumed to be nonzero.
c
c     # NOTE: This routine is obsolete.  Instead of using limiter.f,
c     # which calls philim.f for every wave, it is more efficient to 
c     # use inlinelimiter.f, which eliminates all these function calls
c     # to philim.  If you wish to change the limiter function and are
c     # using inlinelimiter.f, the formulas must be changed in that routine.
c
      r = b/a
      go to (10,20,30,40,50) meth

c
   10 continue
c     --------
c     # minmod
c     --------
      philim = dmax1(0.d0, dmin1(1.d0, r))
      return
c
   20 continue
c     ----------
c     # superbee
c     ----------
      philim = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
      return
c
   30 continue
c     ----------
c     # van Leer
c     ----------
      philim = (r + dabs(r)) / (1.d0 + dabs(r))
      return
c
   40 continue
c     ------------------------------
c     # monotinized centered 
c     ------------------------------
      c = (1.d0 + r)/2.d0
      philim = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
      return
c
   50 continue
c     ------------------------------
c     # Beam-Warming
c     ------------------------------
      philim = r

      return
      end
