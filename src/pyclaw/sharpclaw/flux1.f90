! ===================================================================
subroutine flux1(q1d,dq1d,aux,dt,cfl,t,ixyz,num_aux,num_eqn,mx,num_ghost,maxnx,rp,tfluct)
! ===================================================================
!
!     # Evaluate (delta t) * dq(t)/dt
!
!     The following variables, from the workspace module, are used locally:
!     amdq, apdq, amdq2, apdq2, wave, s 
!
!     amdq(num_eqn,1-num_ghost:mx+num_ghost) = left-going flux-differences
!     apdq(num_eqn,1-num_ghost:mx+num_ghost) = right-going flux-differences
!        e.g. amdq(m,i) = m'th component of A^- \Delta q from i'th Riemann
!                         problem (between cells i-1 and i).
!
!     wave(1-num_ghost:mx+num_ghost, num_eqn, num_waves) = waves from solution of
!                                           Riemann problems,
!            wave(m,mw,i) = mth component of jump in q across
!                           wave in family mw in Riemann problem between
!                           states i-1 and i.
!
!     s(num_waves,1-num_ghost:mx+num_ghost) = wave speeds,
!            s(mw,i) = speed of wave in family mw in Riemann problem between
!                      states i-1 and i.
!
!     Note that mx must be the size of the patch for the dimension corresponding
!     to the value of ixyz.
!
!     t is the time at which we want to evaluate dq/dt, which may not
!      be the current simulation time
! ===================================================================

    USE reconstruct
    USE workspace
    USE ClawParams

    implicit none

    ! Input (dummy) variables
    external :: rp, tfluct
    integer, intent(in) :: num_aux, num_eqn, num_ghost, maxnx, mx, ixyz
    double precision, intent(in) :: q1d(num_eqn,1-num_ghost:mx+num_ghost)
    double precision, intent(inout) :: dq1d(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, target, intent(in) :: aux(num_aux,1-num_ghost:mx+num_ghost)
    double precision, intent(out) :: cfl
    double precision, intent(in) :: t, dt

!f2py intent(in,out) dq1d  
!f2py intent(out) cfl  
!f2py optional dq1d

    ! Local variables
    double precision, pointer :: auxl(:,:), auxr(:,:), qr_shift(:,:), ql_shift(:,:)
    integer :: m, mw, i

! ===================================================================

    dq1d(:,:) = 0.d0

    if (index_capa.gt.0) then
        dtdx = dt / (dx(ixyz)*aux(index_capa,:))
    else
        dtdx = dt / dx(ixyz)
    endif

    select case(lim_type)
        case(1) ! lim_type = 1: 2nd-order TVD reconstruction
        select case(char_decomp)
            case(0)
                ! TVD reconstruction w/o char. decomp.
                call tvd2(q1d,ql,qr,mthlim,num_eqn)
            case(1)
                ! wave-based second order reconstruction
                if (num_dim.eq.1) then
                    call rp(maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,&
                            q1d,q1d,aux,aux,wave,s,amdq,apdq)
                else
                    call rp(ixyz,maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,&
                            q1d,q1d,aux,aux,wave,s,amdq,apdq)
                endif
                ! Need to write a tvd2_fwave routine
                call tvd2_wave(q1d,ql,qr,wave,s,mthlim,num_eqn,num_ghost)
            case(2)
                ! characteristic-wise second order reconstruction
                call evec(mx,num_eqn,num_ghost,mx,q1d,aux,aux,evl,evr)
                call tvd2_char(q1d,ql,qr,mthlim,num_eqn,num_ghost,evl,evr)
        end select
        case(2) ! lim_type = 2: High-order WENO reconstruction
        select case (char_decomp)
            case (0) ! no characteristic decomposition
                call weno_comp(q1d,ql,qr,num_eqn,maxnx,num_ghost)
            case (1) ! wave-based reconstruction
                if (num_dim.eq.1) then
                    call rp(maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,&
                            q1d,q1d,aux,aux,wave,s,amdq,apdq)
                else
                    call rp(ixyz,maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,&
                            q1d,q1d,aux,aux,wave,s,amdq,apdq)
                endif

                if (fwave.eqv. .True.) then
                    call weno5_fwave(q1d,ql,qr,wave,s)
                else
                    call weno5_wave(q1d,ql,qr,wave)
                endif
            case (2) ! characteristic-wise reconstruction
                call evec(mx,num_eqn,num_ghost,mx,q1d,aux,aux,evl,evr)
                call weno5_char(q1d,ql,qr,maxnx,num_eqn,num_ghost,evl,evr)
            case (3) ! transmission-based reconstruction
                call evec(mx,num_eqn,num_ghost,mx,q1d,aux,aux,evl,evr)
                call weno5_trans(q1d,ql,qr,evl,evr)
            case default
                write(*,*) 'ERROR: Unrecognized characteristic decomposition option'
                write(*,*) 'You should set 0<=char_decomp<=3'
                stop
        end select
        case(3)
            call weno5(q1d,ql,qr,num_eqn,maxnx,num_ghost)
        case default
            write(*,*) 'ERROR: Unrecognized limiter type option'
            write(*,*) 'You should set 1<=lim_type<=3'
            stop
      end select


    ! solve Riemann problem at each interface 
    ! -----------------------------------------
    if (num_dim.eq.1) then
        call rp(maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,&
                ql,qr,aux,aux,wave,s,amdq,apdq)
    else
        call rp(ixyz,maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,&
                ql,qr,aux,aux,wave,s,amdq,apdq)
    endif

    ! compute maximum wave speed:
    cfl = 0.d0
    do mw=1,num_waves
        do i=1,mx+1
            ! if s>0 use dtdx(i) to compute CFL,
            ! if s<0 use dtdx(i-1) to compute CFL:
            cfl = dmax1(cfl, dtdx(i)*s(mw,i), -dtdx(i-1)*s(mw,i))
        enddo
    enddo

    ! Find total fluctuation within each cell
    if (tfluct_solver .eqv. .True.) then
        ! tfluct should be a special solver that uses the parameters aux(i)
        ! to solve a Riemann problem with left state ql(i)
        ! and right state qr(i), and returns a total fluctuation in amdq2
        ! NOTE that here amdq2 is really a total fluctuation (should be
        ! called adq); we do it this way just to avoid declaring more storage
        if (num_dim.eq.1) then
            call tfluct(maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,&
                        ql,qr,aux,aux,amdq2)
        else
            call tfluct(ixyz,maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,&
                        ql,qr,aux,aux,amdq2)
        endif

        ! Modify q using fluctuations:
        ! Note this may not correspond to a conservative flux-differencing
        ! for equations not in conservation form. It is conservative if
        ! adq = f(qr(i)) - f(ql(i)).

        forall (i=1:mx, m=1:num_eqn)
            dq1d(m,i) = dq1d(m,i) - dtdx(i)*(apdq(m,i) + &
                            amdq2(m,i) + amdq(m,i+1))
        end forall

    else
        ! Or we can just swap things around and use the usual Riemann solver
        ! This may be more convenient, but is less efficient. 
        
        qr_shift => ql(:,2-num_ghost:mx+num_ghost)
        ql_shift => qr

        auxr => aux(:,2-num_ghost:mx+num_ghost)
        auxl => aux
        
        if (num_dim.eq.1) then
            call rp(maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,&
                    ql_shift,qr_shift,auxl,auxr,wave,s,amdq2,apdq2)
        else
            call rp(ixyz,maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,&
                    ql_shift,qr_shift,auxl,auxr,wave,s,amdq2,apdq2)
        endif

        forall(i=1:mx, m=1:num_eqn)
            dq1d(m,i) = dq1d(m,i)-dtdx(i)*(amdq(m,i+1)+ &
                        apdq(m,i)+amdq2(m,i)+apdq2(m,i))
        end forall
    endif

end subroutine flux1
