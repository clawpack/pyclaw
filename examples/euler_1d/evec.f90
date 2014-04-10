! ==============================================================
subroutine evec(maxnx,num_eqn,num_ghost,mx,q,auxl,auxr,evl,evr)
! ==============================================================
!
!	Calculation of left and right eigenvectors
!
!   On input, q contains the cell average state vector
!             maxnx contains the number of physical points 
!               (without ghost cells)
!             num_ghost is the number of ghost cells
!             num_eqn the number of equations 
!
!   On output, evl(i) and evr(i) contain left/right eigenvectors
!   at interface i-1/2


    implicit double precision (a-h,o-z)
    integer,          intent(in) :: maxnx, num_eqn, num_ghost, mx
    double precision, intent(in) :: q(num_eqn,maxnx+2*num_ghost)
    double precision, intent(out) :: evl(num_eqn,num_eqn,maxnx+2*num_ghost)
    double precision, intent(out) :: evr(num_eqn,num_eqn,maxnx+2*num_ghost)

!   local storage
!   ---------------
    double precision :: u, p, enth, c2, c, det_evr
    common /cparam/ gamma1
    integer :: mx2

    mx2 = size(q,2);

    matrix_type = 2
    ! matrix_type:
    ! 0: computed at cell averages
    ! 1: arithmetic mean of cell averages
    ! 2: Roe mean of cell averages

    do 20 i=2,mx2
        ! Compute velocity, speed and enthalpy:
        select case(matrix_type)
            case(0)
                u = q(2,i) / q(1,i)
                p = gamma1*(q(3,i) - 0.5d0*(q(2,i)**2)/q(1,i))
                enth = (q(3,i)+p) / q(1,i)
                c2 = gamma1*(enth - .5d0*u**2)
                c = dsqrt(c2)
            case(1)
                ul = q(2,i-1) / q(1,i-1)
                ur = q(2,i) / q(1,i)
                u = .5d0*(ul + ur)
                pl = gamma1*(q(3,i-1) - 0.5d0*(q(2,i-1)**2)/q(1,i-1))
                pr = gamma1*(q(3,i) - 0.5d0*(q(2,i)**2)/q(1,i))
                enthl = (q(3,i-1)+pl) / q(1,i-1)
                enthr = (q(3,i)+pr) / q(1,i)
                c2 = .5d0*gamma1*((enthl - .5d0*ul**2) + (enthr - .5d0*ur**2))
                c = dsqrt(c2)
            case(2)
                rhsqrtl = dsqrt(q(1,i-1))
                rhsqrtr = dsqrt(q(1,i  ))
                pl = gamma1*(q(3,i-1) - 0.5d0*(q(2,i-1)**2)/q(1,i-1))
                pr = gamma1*(q(3,i  ) - 0.5d0*(q(2,i  )**2)/q(1,i  ))
                rhsq2 = rhsqrtl + rhsqrtr
                u = (q(2,i-1)/rhsqrtl + q(2,i)/rhsqrtr) / rhsq2
                enth = (((q(3,i-1)+pl)/rhsqrtl + (q(3,i)+pr)/rhsqrtr)) / rhsq2
                c2 = gamma1*(enth - .5d0*u**2)
                c = dsqrt(c2)
        end select


        ! Construct matrix of right eigenvectors
        !      _                     _ 
        !     |                       |
        !     |   1      1       1    |
        !     |                       |
        ! R = |  u-c     u      u+c   |
        !     |                       |
        !     |  H-uc   u^2/2   H+uc  |
        !     |_                     _|

        evr(1,1,i) = 1.d0 
        evr(2,1,i) = u - c
        evr(3,1,i) = enth - u*c

        evr(1,2,i) = 1.d0 
        evr(2,2,i) = u 
        evr(3,2,i) = .5d0*u**2

        evr(1,3,i) = 1.d0 
        evr(2,3,i) = u + c
        evr(3,3,i) = enth + u*c


        ! Construct matrix of left eigenvectors
        !
        ! gamma1 = gamma - 1
        !                          _                                       _ 
        !                         |                                         |
        !                         |  uc/gamma1+u^2/2    -c/gamma1-u     1   |
        !                         |                                         |
        ! R^{-1} =  gamma1/(2c^2) |  2(H-u^2)           2u             -2   |
        !                         |                                         |
        !                         |  -uc/gamma1+u^2/2   c/gamma1-u      1   |
        !                         |_                                       _|

        det_evr_inv = 0.5d0 * gamma1 / c**2

        evl(1,1,i) = (c*u/gamma1 + .5d0*u**2) * det_evr_inv
        evl(2,1,i) = 2.d0*(enth-u**2) * det_evr_inv
        evl(3,1,i) = (-c*u/gamma1 + .5d0*u**2) * det_evr_inv

        evl(1,2,i) = -(c/gamma1 + u) * det_evr_inv
        evl(2,2,i) = 2.d0*u * det_evr_inv
        evl(3,2,i) = (c/gamma1 - u) * det_evr_inv

        evl(1,3,i) = det_evr_inv
        evl(2,3,i) = -2.d0 * det_evr_inv
        evl(3,3,i) = det_evr_inv

    20 enddo

    return
end subroutine evec

