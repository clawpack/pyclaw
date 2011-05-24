module ClawParams

! Problem setup:
  double precision, allocatable :: xlower(:),xupper(:),dx(:)
  integer :: ndim, mwaves, mcapa

! Method-related parameters:
  integer :: char_decomp,lim_type,multid_recon
  integer, allocatable :: mthlim(:)
  logical :: fwave, tfluct_solver

contains

    subroutine alloc_clawparams()

        allocate(xlower(ndim))
        allocate(xupper(ndim))
        allocate(dx(ndim))
        allocate(mthlim(mwaves))

    end subroutine alloc_clawparams

    subroutine dealloc_clawparams()

        deallocate(xlower)
        deallocate(xupper)
        deallocate(dx)
        deallocate(mthlim)

    end subroutine dealloc_clawparams

end module
