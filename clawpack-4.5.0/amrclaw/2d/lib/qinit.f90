! Default qinit file
subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    implicit none
    
    integer, intent(in) :: maxmx,maxmy,meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy
    double precision, intent(inout) :: q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
    double precision, intent(inout) :: aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

end subroutine qinit