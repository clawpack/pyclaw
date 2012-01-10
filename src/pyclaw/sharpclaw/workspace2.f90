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

    subroutine alloc_workspace(maxnx,mbc,meqn,mwaves,char_decomp)
        integer,intent(in) :: maxnx,mbc,meqn,mwaves,char_decomp

        allocate(amdq(meqn,1-mbc:maxnx+mbc))
        allocate(apdq(meqn,1-mbc:maxnx+mbc))
        allocate(amdq2(meqn,1-mbc:maxnx+mbc))
        allocate(apdq2(meqn,1-mbc:maxnx+mbc))
        allocate(ql(meqn,1-mbc:maxnx+mbc))
        allocate(qr(meqn,1-mbc:maxnx+mbc))
        allocate(wave(meqn,mwaves,1-mbc:maxnx+mbc))
        allocate(s(mwaves,1-mbc:maxnx+mbc))
        allocate(dtdx(1-mbc:maxnx+mbc))

        if (char_decomp>1) then
            allocate(evl(meqn,meqn,1-mbc:maxnx+mbc))
            allocate(evr(meqn,meqn,1-mbc:maxnx+mbc))
        endif
            
        work_alloc = .True.

    end subroutine alloc_workspace


    subroutine dealloc_workspace(char_decomp)

        deallocate(amdq)
        deallocate(apdq)
        deallocate(amdq2)
        deallocate(apdq2)
        deallocate(ql)
        deallocate(qr)
        deallocate(wave)
        deallocate(s)
        deallocate(dtdx)

        if (char_decomp>1) then
            deallocate(evl)
            deallocate(evr)
        endif

        work_alloc = .True.

    end subroutine dealloc_workspace

end module
