module rp1_advection
  use iso_c_binding
  implicit none

  type :: advection_context
     real(c_double) :: u
  end type advection_context

contains

  subroutine new(u, context) bind(C, name='rp1_advection_new')
    real(c_double), intent(in)       :: u
    type(c_ptr), intent(out)         :: context

    type(advection_context), pointer :: param

    allocate(param)
    param%u = u
    context = c_loc(param)
  end subroutine new

  subroutine delete(context) bind(C, name='rp1_advection_delete')
    type(c_ptr), intent(in)          :: context

    type(advection_context), pointer :: param

    call c_f_pointer(context, param)
    deallocate(param)
  end subroutine delete

  ! only 1 wave, no need for an mwave parameter
  ! only 1 equation, no need for an meqn parameter

  subroutine advection(mx, mbc, ql, qr, s, wave, amdq, apdq, context) &
       bind(C, name='rp1_advection_c')
    ! C interfacable version of Riemann solver for advection
    !
    !
    ! inputs
    ! * ql  - data at left edge of each cell
    ! * qr  - data at right edge of each cell
    ! * mx  - number of physical grid points
    ! * mbc - number of boundary layer grid points
    ! outputs
    ! * s
    ! * wave
    ! * amdq
    ! * apdq
    ! context
    ! * u - constant advection speed
    !

    ! Solve Riemann problems for the 1D advection equation q_t + u*q_x = 0.
    ! for constant advection velocity u, passed in common block.
    !
    ! The advection speed u is passed in the derived type solver data
    ! On input, ql contains the state vector at the left edge of each cell
    !           qr contains the state vector at the right edge of each cell
    ! On output, wave contains the waves,
    !            s the speeds,
    !            amdq the  left-going flux difference  A^- \Delta q
    !            apdq the right-going flux difference  A^+ \Delta q
    !
    ! Note that the i'th Riemann problem has left state qr(i-1,:) and
    !                                    right state ql(i,:) From the
    !                                    basic clawpack routine step1, rp
    !                                    is called with ql = qr = q.

    integer(c_int), intent(in), value :: mx, mbc
    real(c_double), intent(in) :: &
         ql(1, 1-mbc:mx+mbc), &
         qr(1, 1-mbc:mx+mbc)
    type(c_ptr), intent(in), value :: context

    real(c_double), intent(out) :: &
         s(1, 1-mbc:mx+mbc), &
         wave(1, 1, 1-mbc:mx+mbc), &
         apdq(1, 1-mbc:mx+mbc), &
         amdq(1, 1-mbc:mx+mbc)

    type(advection_context), pointer :: param
    real(c_double) :: u
    integer :: i

    call c_f_pointer(context, param)
    u = param%u

    do i=2-mbc,mx+mbc
       ! Compute the wave
       wave(1,1,i) = ql(1,i) - qr(1,i-1)

       ! Speeds are all the same
       s(1,i) = u

       ! Set fluctuations
       amdq(1,i) = dmin1(u, 0.d0) * wave(1,1,i)
       apdq(1,i) = dmax1(u, 0.d0) * wave(1,1,i)
    end do

  end subroutine advection

end module rp1_advection
