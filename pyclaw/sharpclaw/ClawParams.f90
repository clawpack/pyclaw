module ClawParams
! ===================================================================
! This module holds user options and grid information that are used
! by the SharpClaw solvers in any number of dimensions.
! ===================================================================

    ! Problem setup:
    double precision, allocatable :: xlower(:), xupper(:), dx(:)
    integer :: num_dim, num_waves, index_capa

    ! Method-related parameters:
    integer :: char_decomp, lim_type, weno_order
    integer, allocatable :: mthlim(:)
    logical :: fwave, tfluct_solver

contains

    subroutine alloc_clawparams()

        allocate(xlower(num_dim))
        allocate(xupper(num_dim))
        allocate(dx(num_dim))
        allocate(mthlim(num_waves))

    end subroutine alloc_clawparams

    subroutine dealloc_clawparams()

        if (allocated(xlower)) then
            deallocate(xlower)
        endif

        if (allocated(xupper)) then
            deallocate(xupper)
        endif

        if (allocated(dx)) then
            deallocate(dx)
        endif

        if (allocated(mthlim)) then
            deallocate(mthlim)
        endif

    end subroutine dealloc_clawparams

end module
