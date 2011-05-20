module workspace

  double precision, allocatable :: amdq(:,:),apdq(:,:),dtdx(:)
  double precision, allocatable :: amdq2(:,:),apdq2(:,:)
  double precision, allocatable :: qr(:,:),ql(:,:)
  double precision, allocatable :: s(:,:),wave(:,:,:)
  double precision, allocatable :: evl(:,:,:),evr(:,:,:)
  ! For 2D:
  double precision, allocatable :: aux2(:,:)
  double precision, allocatable :: dtdx1d(:), dtdy1d(:)
  ! For multidimensional reconstruction
  double precision, allocatable :: qgauss(:,:,:,:), q1dgauss(:,:,:)
  double precision, allocatable :: qlgauss(:,:,:), qrgauss(:,:,:)

  logical :: work_alloc = .False.

contains

    subroutine alloc_workspace(mx,my,mbc,meqn,mwaves)
        integer,intent(in) :: mx,my,mbc,meqn,mwaves

        maxnx = max(mx,my)

        allocate(amdq(meqn,1-mbc:maxnx+mbc))
        allocate(apdq(meqn,1-mbc:maxnx+mbc))
        allocate(amdq2(meqn,1-mbc:maxnx+mbc))
        allocate(apdq2(meqn,1-mbc:maxnx+mbc))
        allocate(ql(meqn,1-mbc:maxnx+mbc))
        allocate(qr(meqn,1-mbc:maxnx+mbc))
        allocate(wave(meqn,mwaves,1-mbc:maxnx+mbc))
        allocate(s(mwaves,1-mbc:maxnx+mbc))

        allocate(dtdx(1-mbc:maxnx+mbc))

        work_alloc = .True.

    end subroutine alloc_workspace



end module
