!=====================================================
subroutine evec(maxnx,num_eqn,num_ghost,mx,q,auxl,auxr,evl,evr)
!=====================================================
!
!	# Calculation of left and right eigenvectors
!   # dummy file
!
    implicit double precision (a-h,o-z)
	integer,          intent(in) :: maxnx, num_eqn, num_ghost, mx
	double precision, intent(in) :: q(num_eqn,1-num_ghost:maxnx+num_ghost)
	double precision, intent(out) :: evl(num_eqn,num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(out) :: evr(num_eqn,num_eqn,1-num_ghost:maxnx+num_ghost)

!   # local storage
!   ---------------
	double precision :: u, p, enth, c2, c, det_evl
	double precision :: A(3,3), B(3,3), pro(3,3)
	common /cparam/ gamma1
	integer :: mx2

!   # Construct identity matrices	
    !do 20 m=1-num_ghost,maxnx+num_ghost
    !    do 30 j=1,num_eqn
    !        do 40 i=1,num_eqn
    !            evl(i,j,m) = (i/j)*(j/i)
	!            evr(i,j,m) = (i/j)*(j/i)
    !        40 end do
    !    30 end do
	!20 end do

	do 20 i=1-num_ghost,maxnx+num_ghost
!       # Compute velocity, speed and enthalpy:
        u = q(2,i) / q(1,i)
        p = gamma1*(q(3,i) - 0.5d0*(q(2,i)**2)/q(1,i))
        enth = (q(3,i)+p) / q(1,i)
        c2 = gamma1*(enth - .5d0*u**2)
        c = dsqrt(c2)

!       # Construct matrix of right eigenvectors

        evr(1,1,i) = 1.d0 
        evr(2,1,i) = u - c
        evr(3,1,i) = enth - u*c

        evr(1,2,i) = 1.d0 
        evr(2,2,i) = u 
        evr(3,2,i) = .5d0*u**2

        evr(1,3,i) = 1.d0 
        evr(2,3,i) = u + c
        evr(3,3,i) = enth + u*c

!       # Construct matrix of left eigenvectors

        det_evr = 2.d0*c**3 / gamma1

        evl(1,1,i) = (c*u/gamma1 + .5d0*u**2) * c / det_evr
        evl(1,2,i) = 2.d0*(enth-u**2) * c / det_evr
        evl(1,3,i) = (-c*u/gamma1 + .5d0*u**2) * c / det_evr

        evl(2,1,i) = -(c/gamma1 + u) * c / det_evr
        evl(2,2,i) = 2.d0*u * c / det_evr
        evl(2,3,i) = (c/gamma1 - u) * c / det_evr

        evl(3,1,i) = 1.d0 * c / det_evr
        evl(3,2,i) = -2.d0 * c / det_evr
        evl(3,3,i) = 1.d0 * c / det_evr

!       check that evl*evr = I
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
    20 end do

    return
end subroutine evec

