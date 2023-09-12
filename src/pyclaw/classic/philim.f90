! =====================================================
double precision function philim(a,b,meth)
! =====================================================
    implicit double precision (a-h,o-z)

! Compute a limiter based on wave strengths a and b.
! meth determines what limiter is used.
! a is assumed to be nonzero.

! NOTE: This routine is obsolete.  Instead of using limiter.f,
! which calls philim.f for every wave, it is more efficient to
! use inlinelimiter.f, which eliminates all these function calls
! to philim.  If you wish to change the limiter function and are
! using inlinelimiter.f, the formulas must be changed in that routine.

    r = b/a
    select case (meth)

    case (1) ! minmod
        philim = dmax1(0.d0, dmin1(1.d0, r))

    case (2) ! superbee
        philim = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))

    case (3) ! van Leer
        philim = (r + dabs(r)) / (1.d0 + dabs(r))

    case (4) ! monotonized centered
        c = (1.d0 + r)/2.d0
        philim = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))

    case (5) ! Beam-Warming
        philim = r

    end select
    return
end function
