module workspace
! ===================================================================
! This module holds arrays used as work space to hold intermediate
! values in the SharpClaw solvers.
! ===================================================================

    double precision, allocatable  :: dtdx(:)
    double precision, target, allocatable  :: ql(:,:), qr(:,:)
    double precision, allocatable  :: evl(:,:,:), evr(:,:,:)
    double precision, allocatable  :: amdq(:,:), apdq(:,:), amdq2(:,:), apdq2(:,:)
    double precision, allocatable  :: wave(:,:,:), s(:,:)
    logical :: work_alloc = .False.

contains

    subroutine alloc_workspace(maxnx,num_ghost,num_eqn,num_waves,char_decomp)
        integer,intent(in) :: maxnx,num_ghost,num_eqn,num_waves,char_decomp

        allocate(ql(num_eqn,1-num_ghost:maxnx+num_ghost))
        allocate(qr(num_eqn,1-num_ghost:maxnx+num_ghost))
        allocate(amdq(num_eqn,1-num_ghost:maxnx+num_ghost))
        allocate(apdq(num_eqn,1-num_ghost:maxnx+num_ghost))
        allocate(amdq2(num_eqn,1-num_ghost:maxnx+num_ghost))
        allocate(apdq2(num_eqn,1-num_ghost:maxnx+num_ghost))
        allocate(wave(num_eqn,num_waves,1-num_ghost:maxnx+num_ghost))
        allocate(s(num_waves,1-num_ghost:maxnx+num_ghost))
        allocate(dtdx(1-num_ghost:maxnx+num_ghost))

        if (char_decomp>1) then
            allocate(evl(num_eqn,num_eqn,1-num_ghost:maxnx+num_ghost))
            allocate(evr(num_eqn,num_eqn,1-num_ghost:maxnx+num_ghost))
        endif

        work_alloc = .True.

    end subroutine alloc_workspace

    subroutine dealloc_workspace(char_decomp)

        integer :: char_decomp

        if (allocated(ql)) then
            deallocate(ql)
            deallocate(qr)
            deallocate(amdq)
            deallocate(apdq)
            deallocate(amdq2)
            deallocate(apdq2)
            deallocate(wave)
            deallocate(s)
            deallocate(dtdx)

            if (char_decomp>1) then
                deallocate(evl)
                deallocate(evr)
            endif
        endif

        work_alloc = .False.

    end subroutine dealloc_workspace


end module workspace
