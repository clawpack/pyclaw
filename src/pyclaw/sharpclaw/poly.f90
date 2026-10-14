! ===================================================================
module poly
! ===================================================================
!   Linear (unlimited) polynomial reconstruction of interface values
!   from cell averages, for the SharpClaw finite-volume solvers.
!
!   q(:,i) holds the average of q over cell i. For each cell, ql(:,i) is
!   the reconstructed value at the left face x_(i-1/2) and qr(:,i) the
!   value at the right face x_(i+1/2). The weights are exact: they are
!   obtained by differentiating the polynomial that interpolates the
!   primitive function of q at the cell faces, so they are correct for
!   cell averages (not point values). They are applied component-wise.
!
!   lim_type selects the stencil; order = solver.reconstruction_order,
!   an even number N in {4, 6, 8, 10}. With k = N/2:
!
!     4 ('poly')               ql(i) uses cells i-k..i+k-1 (N cells,
!                              centred on the face). Order N. ql(i+1)
!                              equals qr(i), so the scheme adds no
!                              numerical dissipation.
!     5 ('poly_dissipative')   ql(i) uses cells i-k+1..i+k (the central
!                              stencil shifted one cell upwind for each
!                              face state). Order N, dissipative.
!     6 ('poly_cell_centered') ql(i) uses cells i-k..i+k (N+1 cells,
!                              centred on the cell, polynomial degree N).
!                              Order N+1, dissipative; equal to WENO of
!                              order N+1 with its ideal linear weights.
!
!   qr uses the mirror image of the ql stencil. Every variant reaches k
!   cells from the cell being reconstructed, so num_ghost = k + 1.
!
!   Coefficients generated with exact rational arithmetic.
! ===================================================================

    implicit none

contains

    ! ===================================================================
    subroutine poly_weights(lim_type, order, lo, hi, w, ok)
    ! ===================================================================
    !   Weights for the left-face value ql(i) = sum_j w(j) * q(i+j),
    !   j = lo..hi.

        integer, intent(in) :: lim_type, order
        integer, intent(out) :: lo, hi
        double precision, intent(out) :: w(-5:5)
        logical, intent(out) :: ok

        ok = .true.
        lo = 0; hi = 0
        w = 0.d0

        select case (lim_type)
        case (4) ! central, non-dissipative
            select case (order)
            case (4)
                lo = -2; hi = 1
                w(-2:1) = [-1.d0, 7.d0, 7.d0, -1.d0] / 12.d0
            case (6)
                lo = -3; hi = 2
                w(-3:2) = [1.d0, -8.d0, 37.d0, 37.d0, -8.d0, 1.d0] / 60.d0
            case (8)
                lo = -4; hi = 3
                w(-4:3) = [-3.d0, 29.d0, -139.d0, 533.d0, 533.d0, -139.d0, 29.d0, -3.d0] / 840.d0
            case (10)
                lo = -5; hi = 4
                w(-5:4) = [2.d0, -23.d0, 127.d0, -473.d0, 1627.d0, &
                      1627.d0, -473.d0, 127.d0, -23.d0, 2.d0] / 2520.d0
            case default
                ok = .false.
            end select
        case (5) ! upwind-biased, dissipative
            select case (order)
            case (4)
                lo = -1; hi = 2
                w(-1:2) = [3.d0, 13.d0, -5.d0, 1.d0] / 12.d0
            case (6)
                lo = -2; hi = 3
                w(-2:3) = [-2.d0, 22.d0, 57.d0, -23.d0, 7.d0, -1.d0] / 60.d0
            case (8)
                lo = -3; hi = 4
                w(-3:4) = [5.d0, -55.d0, 365.d0, 743.d0, -307.d0, 113.d0, -27.d0, 3.d0] / 840.d0
            case (10)
                lo = -4; hi = 5
                w(-4:5) = [-3.d0, 37.d0, -233.d0, 1207.d0, 2131.d0, &
                      -893.d0, 367.d0, -113.d0, 22.d0, -2.d0] / 2520.d0
            case default
                ok = .false.
            end select
        case (6) ! cell-centred, dissipative
            select case (order)
            case (4)
                lo = -2; hi = 2
                w(-2:2) = [-3.d0, 27.d0, 47.d0, -13.d0, 2.d0] / 60.d0
            case (6)
                lo = -3; hi = 3
                w(-3:3) = [4.d0, -38.d0, 214.d0, 319.d0, -101.d0, 25.d0, -3.d0] / 420.d0
            case (8)
                lo = -4; hi = 4
                w(-4:4) = [-5.d0, 55.d0, -305.d0, 1375.d0, &
                      1879.d0, -641.d0, 199.d0, -41.d0, 4.d0] / 2520.d0
            case (10)
                lo = -5; hi = 5
                w(-5:5) = [12.d0, -153.d0, 947.d0, -4003.d0, 15797.d0, &
                      20417.d0, -7303.d0, 2597.d0, -703.d0, 122.d0, -10.d0] / 27720.d0
            case default
                ok = .false.
            end select
        case default
            ok = .false.
        end select

    end subroutine poly_weights

    ! ===================================================================
    subroutine poly_recon(q,ql,qr,num_eqn,maxnx,num_ghost,lim_type,order)
    ! ===================================================================

        integer,          intent(in) :: num_eqn, maxnx, num_ghost, lim_type, order
        double precision, intent(in) :: q(num_eqn,maxnx+2*num_ghost)
        double precision, intent(out) :: ql(num_eqn,maxnx+2*num_ghost),qr(num_eqn,maxnx+2*num_ghost)

        double precision :: w(-5:5)
        integer :: lo, hi, reach, i, j, mx2
        logical :: ok

        call poly_weights(lim_type, order, lo, hi, w, ok)
        if (.not. ok) then
            print *, 'ERROR: polynomial reconstruction needs lim_type 4, 5 or 6'
            print *, '       and an even reconstruction_order between 4 and 10.'
            stop
        endif

        mx2 = size(q,2)
        reach = max(-lo, hi)
        if (reach >= num_ghost) then
            print *, 'ERROR: polynomial reconstruction needs num_ghost >= ', reach+1
            stop
        endif

        ql = 0.d0
        qr = 0.d0
        do i = 1+reach, mx2-reach
            do j = lo, hi
                ql(:,i) = ql(:,i) + w(j)*q(:,i+j)
                ! qr mirrors the ql stencil about cell i
                qr(:,i) = qr(:,i) + w(j)*q(:,i-j)
            end do
        end do

    end subroutine poly_recon

end module poly
