! Default qinit function
subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)

    implicit none
    integer, intent(in) :: maxmx,meqn,mbc,mx,maux
    double precision, intent(in) :: xlower,dx
    double precision, intent(inout) :: q(1-mbc:maxmx+mbc,meqn)
    double precision, intent(inout) :: aux(1-mbc:maxmx+mbc,maux)

end subroutine qinit