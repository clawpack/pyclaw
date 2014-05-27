! ============================================================================
subroutine tfluct(ixyz,maxnx,num_eqn,num_waves,num_ghost,mx,ql,qr,auxl,auxr,adq)
! ============================================================================
!   
!   "Internal" Riemann solver for the euler equations in 1D.
!   Te riemann problem is solved by assuming a discontinuity at the
!   center of the i'th cell.
!
!   On input, q contains the cell average state vector
!             ql contains the state vector at the left edge of each cell
!             qr contains the state vector at the right edge of each cell
!             auxl contains the auxiliary vector at the left edge of each cell
!             auxr contains the state vector at the right edge of each cell
!             maxnx contains the number of physical points 
!               (without ghost cells)
!             num_ghost is the number of ghost cells
!             num_eqn the number of equations 
!
!   On output, adq contains the decomposition of the flux difference
!              f(qr(i)) - f(ql(i)).
!
!   For the Euler equations q = (q1, q2, q3)^T, 
!   where q1 = rho, q2 = rho*u, q3 = E and
!              _                                              _  
!             |                       q2                       |
!   f(q) =    |        0.5(3-gamma)q2^2/q1 + (gamma-1)q3       |
!             |_  ( gamma*q3 - 0.5(gamma-1)q2^2/q1 ) * q2/g1  _|


    implicit double precision (a-h,o-z)
    integer,          intent(in) :: maxnx, num_eqn, num_waves, num_ghost, mx
    double precision, intent(in) :: ql(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in) :: qr(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(out) :: adq(num_eqn,1-num_ghost:maxnx+num_ghost)

!   local storage
!   ---------------
    common /cparam/  gamma1
    double precision :: gamma

    gamma = gamma1 + 1.d0

    do i = 1,mx
        adq(1,i) = qr(2,i) - ql(2,i)
        adq(2,i) = 0.5d0*(3.d0-gamma)*(qr(2,i)**2/qr(1,i) - ql(2,i)**2/ql(1,i)) &
                    + gamma1*(qr(3,i) - ql(3,i))
        adq(3,i) = (gamma*qr(3,i)-0.5d0*gamma1*qr(2,i)**2/qr(1,i))*qr(2,i)/qr(1,i) & 
                    - (gamma*ql(3,i)-0.5d0*gamma1*ql(2,i)**2/ql(1,i))*ql(2,i)/ql(1,i)
    enddo

    return
end subroutine tfluct
