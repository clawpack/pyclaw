! ============================================================================
!  Program:     AMRClaw
!  File:        resize_storage.f90
!  Created:     2009-01-21
!  Author:      Kyle Mandli and Marsha Berger
! ============================================================================
!  Description:  Resize the alloc array for AMR storage
! ============================================================================


! NOTE:  Older f90 compilers (e.g. gfortran prior to 4.2?)
! may not implement move_alloc.  If this fails, you may need to use
! resize_storage_static.f90 instead of this routine and set the
! allocation large enough in init_alloc.f90 to avoid running out of space.


subroutine resize_storage(new_size)
    
    use mem_storage
    implicit none
    
    integer, intent(in) :: new_size
    double precision, pointer, dimension(:) :: alloc
    integer :: memsize, status
    common /calloc/ alloc,memsize
    
    double precision, allocatable, target, dimension(:) :: new_storage
    

    if (memsize < new_size) then
        print *, "Expanding storage from ", memsize," to ", new_size
        allocate(new_storage(new_size),STAT=status)
        if (status > 0) then
            print *, 'Allocation failed, status = ', status
            stop
        endif
        new_storage(1:memsize) = storage

        call move_alloc(new_storage,storage)

        alloc => storage
        memsize = new_size
    else
        print *,'new_size < memsize,'
        stop
    endif
    
end subroutine resize_storage
