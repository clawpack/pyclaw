module ClawParams

! Dimensions of arrays:
  integer :: ndim,mcapa,maux,meqn,mwaves,mbc,maxnx
  integer, allocatable :: nx(:)

! Problem setup:
  double precision, allocatable :: xlower(:),xupper(:),dx(:)
  integer, allocatable :: mthbc(:)

! Method-related parameters:
  integer :: src_term,tfluct_solver,char_decomp,lim_type,multid_recon
  integer, allocatable :: mthlim(:)
  integer :: rord

end module
