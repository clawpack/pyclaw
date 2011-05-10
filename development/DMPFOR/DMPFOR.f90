subroutine DMPFOR(a, nx, ny, nz)

    implicit none
    integer, intent(in) :: nx, ny, nz
    double precision, intent(in) :: a(1:nx,1:ny ,1:nz)
    
    integer :: i,j,k
    print *, "Dump from inside fortran routine. Following fortran memory layout"
    do k = 1,nz
        print *, "k=", k
        do j = 1,ny
            print *, "    j=", j
            
            print *, a(:,j,k)
        
        
        enddo
    enddo

end subroutine DMPFOR

program dmpfor_test
    print *, "inside main"

end program dmpfor_test
