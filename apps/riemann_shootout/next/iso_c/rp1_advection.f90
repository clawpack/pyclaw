use iso_c_binding, only: c_ptr
implicit none

integer, parameter:: dp=kind(0.d0)

! =========================================================
subroutine c_advection_1(mx, mbc, ql, qr, s, wave, amdq, apdq, u)
! =========================================================
! C interfacable version of Riemann solver for advection
!
! THIS CODE IS CURRENTLY BROKEN!
!
! static_data
! * u - constant advection speed
!
! state_in
! * ql  - data at left edge of each cell
! * qr  - data at right edge of each cell
! * mx  - number of physical grid points
! * mbc - number of boundary layer grid points
! state_out
! * s
! * wave
! * amdq
! * apdq
!
! Solve Riemann problems for the 1D advection equation q_t + u*q_x = 0.
! for constant advection velocity u, passed in common block.
!
! The advection speed u is passed in the common block cparam
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

  do i=2-mbc,mx+mbc
     !        # Compute the wave and speed
     wave(i) = ql(i) - qr(i-1)
     s(i) = u
     amdq(i) = dmin1(u, 0.d0) * wave(i)
     apdq(i) = dmax1(u, 0.d0) * wave(i)
  end do

  return
end subroutine advection_1

subroutine advection_1(static, state_in, state_out), bind (c)
  type f_state_in
     integer  :: mx
     integer  :: mbc
     real(dp) :: ql(1-mbc:mx+mbc)
     real(dp) :: qr(1-mbc:mx+mbc)




  real(c_double), intent(in)  :: ql(1-mbc:mx+mbc)
  real(c_double), intent(in)  :: qr(1-mbc:mx+mbc)
  real(c_double), intent(out) :: amdq(1-mbc:mx+mbc)
  real(c_double), intent(out) :: apdq(1-mbc:mx+mbc)
  real(c_double), intent(out) :: s(1-mbc:mx+mbc)
  real(c_double), intent(out) :: wave(1-mbc:mx+mbc)
