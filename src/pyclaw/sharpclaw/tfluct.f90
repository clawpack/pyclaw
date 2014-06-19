! =========================================================
subroutine tfluct(ixyz,maxnx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,amdq2)
! =========================================================
!
!	dummy tfluct function (see example euler_1d for more information)

    implicit double precision (a-h,o-z)

    write(*,*) "Error: method(8)=1, but you have not defined"
    write(*,*) "a function tfluct1()."

    stop

    return
end subroutine tfluct
