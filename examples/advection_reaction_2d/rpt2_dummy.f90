! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
    implicit double precision (a-h,o-z)

!     # Dummy transverse Riemann solver, for use in dimensionally-split algorithm.

    write(*,*) 'Error: Dummy transverse Riemann solver called!'
    return
    end subroutine rpt2
