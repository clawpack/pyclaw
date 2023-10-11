! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================
! Riemann solver for the non-conservative transport equation
!    q_t  +  u(x,y,t)*q_x + v(x,y,t)*q_y = 0
! where u and v are a given velocity field and q has multiple components.
! The velocities are specified at cell edges, so that 

!       u(i,j) = u(x_(i-1/2), y_j)
!       v(i,j) = v(x_i, y_(j-1/2))

! The size of q should be given by meqn.

! waves: 1
! equations: meqn
! aux fields: 2

! Conserved quantities: meqn

! Auxiliary variables:
!         1  x_velocity (u)
!         2  y_velocity (v)

! solve Riemann problems along one slice of data.
! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! On output, wave contains the waves, s the speeds,
! and amdq, apdq the left-going and right-going flux differences,
! respectively.  Note that in this advective form, the sum of
! amdq and apdq is not equal to a difference of fluxes except in the
! case of constant velocities.

! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr

    implicit real*8(a-h,o-z)

    dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
    dimension    s(mwaves, 1-mbc:maxm+mbc)
    dimension   ql(meqn, 1-mbc:maxm+mbc)
    dimension   qr(meqn, 1-mbc:maxm+mbc)
    dimension apdq(meqn,1-mbc:maxm+mbc)
    dimension amdq(meqn,1-mbc:maxm+mbc)
    dimension auxl(maux,1-mbc:maxm+mbc)
    dimension auxr(maux,1-mbc:maxm+mbc)


    do i = 2-mbc, mx+mbc
        wave(:,1,i) = ql(:,i) - qr(:,i-1)
        s(1,i) = auxl(ixy,i)
        ! The flux difference df = s*wave  all goes in the downwind direction:
        amdq(:,i) = dmin1(auxl(ixy,i), 0.d0) * wave(:,1,i)
        apdq(:,i) = dmax1(auxl(ixy,i), 0.d0) * wave(:,1,i)
    end do

    return
    end subroutine rpn2
