! =============================================================================
	subroutine rp1(maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =============================================================================
!
! Riemann problems for the 1D Burgers' equation with entropy fix for 
! transonic rarefaction. See "Finite Volume Method for Hyperbolic Problems",
! R. J. LeVeque.

    
    implicit none 
	
	integer :: maxmx, meqn, mwaves, mbc, mx
	
    double precision :: ql(meqn,1-mbc:maxmx+mbc)
    double precision :: qr(meqn,1-mbc:maxmx+mbc)
    double precision :: auxl(1,1-mbc:maxmx+mbc)
    double precision :: auxr(1,1-mbc:maxmx+mbc)
    double precision :: s(mwaves, 1-mbc:maxmx+mbc)
    double precision :: wave(meqn, mwaves, 1-mbc:maxmx+mbc)
    double precision :: amdq(meqn, 1-mbc:maxmx+mbc)
    double precision :: apdq(meqn, 1-mbc:maxmx+mbc)
    
    integer :: i

    do i=2-mbc,mx+mbc
		wave(1,1,i) = ql(1,i) - qr(1,i-1)
        s(1,i) = 0.5d0 * (qr(1,i-1) + ql(1,i))

        amdq(i,1) = dmin1(s(1,i), 0.d0) * wave(1,1,i)
        apdq(i,1) = dmax1(s(1,i), 0.d0) * wave(1,1,i)

        if (ql(1,i).gt.0.d0 .and. qr(1,i-1).lt.0.d0) then
        	amdq(1,i) = - 1.d0/2.d0 * qr(1,i-1)**2
            apdq(1,i) =   1.d0/2.d0 * ql(1,i)**2
        endif
	enddo

    return
    end subroutine