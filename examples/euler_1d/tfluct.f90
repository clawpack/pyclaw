! ============================================================================
subroutine tfluct1(maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,ql,qr,auxl,auxr,amdq2)
! ============================================================================
!   "Internal" Riemann solver for the euler equations in 1D.
!   The riemann problem is solved by assuming a discontinuity at the
!   center of the i'th cell.
!
!   On input, ql contains the state vector at the left edge of each cell
!             qr contains the state vector at the right edge of each cell
!             auxl contains the auxiliary vector at the left edge of each cell
!             auxr contains the state vector at the right edge of each cell
!             maxnx is the number of physical points (without ghost cells)
!             num_ghost is the number of ghost cells
!             num_eqn is the number of equations
!             ixyz is the dimension index
!             mx is the size of the patch for the dimension corresponding
!               to the value of ixyz
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


    implicit none
    integer,          intent(in)  :: maxnx, mx, num_eqn, num_waves, num_aux, num_ghost
    double precision, intent(in)  :: auxl(num_aux,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  :: auxr(num_aux,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  ::   ql(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  ::   qr(num_eqn,1-num_ghost:maxnx+num_ghost)

    double precision, intent(out) :: amdq2(num_eqn,1-num_ghost:maxnx+num_ghost)

    integer :: i
    double precision :: gamma, gamma1

!   local storage
!   ---------------
    common /cparam/  gamma1
    

    gamma = gamma1 + 1.d0

    do i = 1,mx
        amdq2(1,i) = qr(2,i) - ql(2,i)
        amdq2(2,i) = 0.5d0*(3.d0-gamma)*(qr(2,i)**2/qr(1,i) - ql(2,i)**2/ql(1,i)) &
                    + gamma1*(qr(3,i) - ql(3,i))
        amdq2(3,i) = (gamma*qr(3,i)-0.5d0*gamma1*qr(2,i)**2/qr(1,i))*qr(2,i)/qr(1,i) & 
                    - (gamma*ql(3,i)-0.5d0*gamma1*ql(2,i)**2/ql(1,i))*ql(2,i)/ql(1,i)
    enddo

    return
end subroutine tfluct1

