

!     =====================================================
    double precision function philim(a,b,meth)
!     =====================================================
    implicit double precision (a-h,o-z)

!     # Compute a limiter based on wave strengths a and b.
!     # meth determines what limiter is used.
!     # a is assumed to be nonzero.

!     # NOTE: This routine is obsolete.  Instead of using limiter.f,
!     # which calls philim.f for every wave, it is more efficient to
!     # use inlinelimiter.f, which eliminates all these function calls
!     # to philim.  If you wish to change the limiter function and are
!     # using inlinelimiter.f, the formulas must be changed in that routine.

    r = b/a
    go to (10,20,30,40,50) meth


    10 continue
!     --------
!     # minmod
!     --------
    philim = dmax1(0.d0, dmin1(1.d0, r))
    return

    20 continue
!     ----------
!     # superbee
!     ----------
    philim = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
    return

    30 continue
!     ----------
!     # van Leer
!     ----------
    philim = (r + dabs(r)) / (1.d0 + dabs(r))
    return

    40 continue
!     ------------------------------
!     # monotinized centered
!     ------------------------------
    c = (1.d0 + r)/2.d0
    philim = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
    return

    50 continue
!     ------------------------------
!     # Beam-Warming
!     ------------------------------
    philim = r

    return
    END FUNCTION
