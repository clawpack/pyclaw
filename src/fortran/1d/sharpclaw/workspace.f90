module workspace
    double precision, allocatable  :: ql(:,:), qr(:,:), dtdx(:)
    double precision, allocatable  :: amdq(:,:), apdq(:,:), amdq2(:,:), apdq2(:,:)
    double precision, allocatable  :: wave(:,:,:), s(:,:)
    logical :: work_alloc = .False.

contains

    subroutine alloc_workspace(maxnx,mbc,meqn,mwaves)
        integer,intent(in) :: maxnx,mbc,meqn,mwaves

        allocate(ql(meqn,1-mbc:maxnx+mbc))
        allocate(qr(meqn,1-mbc:maxnx+mbc))
        allocate(amdq(meqn,1-mbc:maxnx+mbc))
        allocate(apdq(meqn,1-mbc:maxnx+mbc))
        allocate(amdq2(meqn,1-mbc:maxnx+mbc))
        allocate(apdq2(meqn,1-mbc:maxnx+mbc))
        allocate(wave(meqn,mwaves,1-mbc:maxnx+mbc))
        allocate(s(mwaves,1-mbc:maxnx+mbc))

        allocate(dtdx(1-mbc:maxnx+mbc))

        work_alloc = .True.

    end subroutine alloc_workspace

    subroutine dealloc_workspace()

        deallocate(ql)
        deallocate(qr)
        deallocate(amdq)
        deallocate(apdq)
        deallocate(amdq2)
        deallocate(apdq2)
        deallocate(wave)
        deallocate(s)

        deallocate(dtdx)

        work_alloc = .False.

    end subroutine dealloc_workspace


end module workspace
