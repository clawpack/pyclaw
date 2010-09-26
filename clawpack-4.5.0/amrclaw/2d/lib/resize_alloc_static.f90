
! ============================================================================
!  Program:     AMRClaw
!  File:        resize_storage_static.f90
!  Created:     2009-10-20
!  Author:      Marsha Berger and Randy LeVeque
! ============================================================================
!  Description: For use with compilers that don't support move_alloc.
!  Halts with an error message.
! ============================================================================

! NOTE:  Older f90 compilers (e.g. gfortran prior to 4.2?)
! may not implement move_alloc and you will need to use this routine
! instead of resize_storage.f90 and set the
! allocation large enough in init_alloc.f90 to avoid running out of space.


subroutine resize_storage(new_size)
    
    use mem_storage
    implicit none
    
    integer, intent(in) :: new_size
    double precision, pointer, dimension(:) :: alloc
    integer :: memsize, status
    logical :: no_move_alloc
    common /calloc/ alloc,memsize

    write(6,*) '*** Ran out of storage for AMR.  '
    write(6,*) '*** Increase memsize in init_alloc_static.f'
    write(6,*) '*** or switch to dynamic memory using init_alloc.f90'
    write(6,*) '*** Current memsize = ',memsize
    write(6,*) '*** Requested new_size = ',new_size
    stop
    
end subroutine resize_storage
