! ============================================================================
!  Program:     AMRClaw
!  File:        restrt_alloc.f90
!  Created:     2009-10-22
!  Author:      Marsha Berger and Randy LeVeque
! ============================================================================
!  Description:  Initialization of alloc storage for restart
! ============================================================================


subroutine restrt_alloc(isize)

    use mem_storage
    implicit none
    integer :: isize

    double precision, pointer, dimension(:) :: alloc
    integer :: memsize
    common  /calloc/ alloc, memsize

    if (.not.allocated(storage)) then
        write(6,*)"allocating ",isize," -sized alloc array"
	memsize = isize
        allocate(storage(memsize))
        alloc => storage
        print *, "Storage allocated of size ",memsize," at restart"
    else
        print *, "Storage already allocated!"
    endif

end subroutine restrt_alloc


