!=====================================================
subroutine evec(maxnx,num_eqn,num_ghost,mx,q,ql,qr,evl,evr,flag)
!=====================================================
!
!	Calculation of left and right eigenvectors
!
    implicit double precision (a-h,o-z)
	integer,          intent(in) :: maxnx, num_eqn, num_ghost, mx
	double precision, intent(in) :: q(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in) :: ql(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in) :: qr(num_eqn,1-num_ghost:maxnx+num_ghost)
	double precision, intent(out) :: evl(num_eqn,num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(out) :: evr(num_eqn,num_eqn,1-num_ghost:maxnx+num_ghost)
    logical, intent(out) :: flag

!   local storage
!   ---------------
	double precision :: u, p, enth, c2, c, det_evl
	double precision :: A(3,3), B(3,3), pro(3,3)
	common /cparam/ gamma1
	integer :: mx2
    
    ! Construct identity matrices	
    !do 20 m=1-num_ghost,maxnx+num_ghost
    !    do 30 j=1,num_eqn
    !        do 40 i=1,num_eqn
    !            evl(i,j,m) = (i/j)*(j/i)
	!            evr(i,j,m) = (i/j)*(j/i)
    !        40 end do
    !    30 end do
	!20 end do


   ! if ((ALL(ql.eq.0d0)) .AND. (ALL(qr.eq.0d0))) then
   !     matrix_type = 0
   !     flag = .False.
   ! else
   !     matrix_type = 3
   !     flag = .True.
   ! endif
	
    ! matrix_type:
    ! 0: computed at cell avarages
    ! 1: arithmetic mean of cell averages
    ! 2: Roe mean of cell averages
    ! 3: Roe averages of ql and qr

    matrix_type = 0

    do 20 i=1-num_ghost,maxnx+num_ghost
        ! Compute velocity, speed and enthalpy:
        select case(matrix_type)
            case(0)
                u = q(2,i) / q(1,i)
                p = gamma1*(q(3,i) - 0.5d0*(q(2,i)**2)/q(1,i))
                enth = (q(3,i)+p) / q(1,i)
                c2 = gamma1*(enth - .5d0*u**2)
                c = dsqrt(c2)
            case(1)
                u = .5d0*(q(2,i-1) / q(1,i-1) + q(2,i) / q(1,i))
                pl = gamma1*(q(3,i-1) - 0.5d0*(q(2,i-1)**2)/q(1,i-1))
                pr = gamma1*(q(3,i) - 0.5d0*(q(2,i)**2)/q(1,i))
                enth = .5d0*((q(3,i-1)+pl) / q(1,i-1) + (q(3,i)+pr) / q(1,i))
                c2 = gamma1*(enth - .5d0*u**2)
                c = dsqrt(c2)
            case(2)
                rhsqrtl = dsqrt(q(1,i-1))
                rhsqrtr = dsqrt(q(1,i))
                pl = gamma1*(q(3,i-1) - 0.5d0*(q(2,i-1)**2)/q(1,i-1))
                pr = gamma1*(q(3,i) - 0.5d0*(q(2,i)**2)/q(1,i))
                rhsq2 = rhsqrtl + rhsqrtr
                u = (q(2,i-1)/rhsqrtl + q(2,i)/rhsqrtr) / rhsq2
                enth = (((q(3,i-1)+pl)/rhsqrtl + (q(3,i)+pr)/rhsqrtr)) / rhsq2
                c2 = gamma1*(enth - .5d0*u**2)
                c = dsqrt(c2)
            case(3)
                rhsqrtl = dsqrt(qr(1,i-1))
                rhsqrtr = dsqrt(ql(1,i))
                pl = gamma1*(qr(3,i-1) - 0.5d0*(qr(2,i-1)**2)/qr(1,i-1))
                pr = gamma1*(ql(3,i) - 0.5d0*(ql(2,i)**2)/ql(1,i))
                rhsq2 = rhsqrtl + rhsqrtr
                u = (qr(2,i-1)/rhsqrtl + ql(2,i)/rhsqrtr) / rhsq2
                enth = (((qr(3,i-1)+pl)/rhsqrtl &
                        + (ql(3,i)+pr)/rhsqrtr)) / rhsq2
                c2 = gamma1*(enth - .5d0*u**2)
                c = dsqrt(c2)
        end select


        ! Construct matrix of right eigenvectors

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

        det_evr = 2.d0*c**3 / gamma1

       ! evl(1,1,i) = (c*u/gamma1 + .5d0*u**2) * c / det_evr
       ! evl(1,2,i) = 2.d0*(enth-u**2) * c / det_evr
       ! evl(1,3,i) = (-c*u/gamma1 + .5d0*u**2) * c / det_evr

       ! evl(2,1,i) = -(c/gamma1 + u) * c / det_evr
       ! evl(2,2,i) = 2.d0*u * c / det_evr
       ! evl(2,3,i) = (c/gamma1 - u) * c / det_evr

       ! evl(3,1,i) = 1.d0 * c / det_evr
       ! evl(3,2,i) = -2.d0 * c / det_evr
       ! evl(3,3,i) = 1.d0 * c / det_evr

        evl(1,1,i) = (c*u/gamma1 + .5d0*u**2) * c / det_evr
        evl(2,1,i) = 2.d0*(enth-u**2) * c / det_evr
        evl(3,1,i) = (-c*u/gamma1 + .5d0*u**2) * c / det_evr

        evl(1,2,i) = -(c/gamma1 + u) * c / det_evr
        evl(2,2,i) = 2.d0*u * c / det_evr
        evl(3,2,i) = (c/gamma1 - u) * c / det_evr

        evl(1,3,i) = 1.d0 * c / det_evr
        evl(2,3,i) = -2.d0 * c / det_evr
        evl(3,3,i) = 1.d0 * c / det_evr

       ! check that evl*evr = I
       ! do m=1,3
       !     do n=1,3
       !          A(m,n) = evl(m,n,i)
       !          B(m,n) = evr(m,n,i)
       !      end do
       ! end do

       ! do m = 1,3
       !     do n = 1,3
       !     summ = 0.d0
       !         do k = 1,3
       !             summ = summ + A(m,k)*B(k,n)
       !         end do
       !     pro(m,n) = summ
       !     end do
       ! end do
       ! 
       ! print *, pro

       ! write( *, * ) 'Press Enter to continue'
       ! read( *, * )
    20 enddo

    return
end subroutine evec

