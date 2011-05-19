module ClawParams

! Problem setup:
  double precision, allocatable :: xlower(:),xupper(:),dx(:)

! Method-related parameters:
  integer :: tfluct_solver,char_decomp,lim_type,multid_recon
  integer, allocatable :: mthlim(:)

contains

    subroutine alloc_clawparams(ndim,mwaves)
        integer, intent(in) :: ndim, mwaves

        allocate(xlower(ndim))
        allocate(xupper(ndim))
        allocate(dx(ndim))
        allocate(mthlim(mwaves))

    end subroutine alloc_clawparams

end module
