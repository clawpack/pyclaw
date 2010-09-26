! ============================================================================
!  Program:     AMRClaw
!  File:        init_alloc.f90
!  Created:     2009-01-21
!  Author:      Kyle Mandli and Marsha Berger
! ============================================================================
!  Description:  Initialization of alloc storage
! ============================================================================

module mem_storage
    double precision, allocatable, target, dimension(:) :: storage
end module mem_storage

subroutine init_alloc()
    
    use mem_storage
    implicit none
    
    double precision, pointer, dimension(:) :: alloc
    integer :: memsize
    common  /calloc/ alloc, memsize
    
    if (.not.allocated(storage)) then
        memsize = 1000000
        allocate(storage(memsize))
        alloc => storage
        print *, "Storage allocated..."
    else
        print *, "Storage already allocated!"
    endif
    
end subroutine init_alloc

