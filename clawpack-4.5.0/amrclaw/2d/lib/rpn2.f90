! Default rpn2 routine
subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)

    implicit none
    
    integer, intent(in) :: ixy,maxm,meqn,mwaves,mbc,mx
    double precision, intent(inout) :: ql(1-mbc:maxm+mbc,meqn)
    double precision, intent(inout) :: qr(1-mbc:maxm+mbc,meqn)
    double precision, intent(inout) :: auxl(1-mbc:maxm+mbc,meqn)
    double precision, intent(inout) :: auxr(1-mbc:maxm+mbc,meqn)
    double precision, intent(out) :: wave(1-mbc:maxm+mbc,meqn,mwaves)
    double precision, intent(out) :: s(1-mbc:maxm+mbc, mwaves)
    double precision, intent(out) :: apdq(1-mbc:maxm+mbc, meqn)
    double precision, intent(out) :: amdq(1-mbc:maxm+mbc, meqn)
    
    wave = 0.d0
    s = 0.d0
    apdq = 0.d0
    amdq = 0.d0
    
end subroutine rpn2
