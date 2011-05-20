module ClawParams

! Problem setup:
  double precision, allocatable :: xlower(:),xupper(:),dx(:)
  integer :: ndim, mwaves, mcapa

! Method-related parameters:
  integer :: tfluct_solver,char_decomp,lim_type,multid_recon
  integer, allocatable :: mthlim(:)

contains

    subroutine alloc_clawparams()

        allocate(xlower(ndim))
        allocate(xupper(ndim))
        allocate(dx(ndim))
        allocate(mthlim(mwaves))

    end subroutine alloc_clawparams

end module
