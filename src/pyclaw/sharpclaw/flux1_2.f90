! ===================================================================
subroutine flux1(q1d,dq1d,aux,dt,cfl,t,ixy,num_aux,num_eqn,mx,num_ghost,maxnx,rpn2)
! ===================================================================
!
!     # Evaluate (delta t) * dq(t)/dt
!
!     SharpClaw
!     Author: David Ketcheson
!
!     amdq, apdq, amdq2, apdq2, wave, and s are used locally:
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
!     to the value of ixy.
!
!     t is the time at which we want to evaluate dq/dt, which may not
!      be the current simulation time
! ===================================================================
!
! Modified: April 26, 2011
! Authors:  David Ketcheson
!           Matteo Parsani
!
! ===================================================================

    USE reconstruct
    USE workspace
    USE ClawParams

    implicit double precision (a-h,o-z)

    integer :: num_aux,num_eqn,num_ghost,maxnx,mx
    double precision :: q1d(num_eqn,1-num_ghost:mx+num_ghost)
    double precision :: dq1d(num_eqn,1-num_ghost:maxnx+num_ghost)
    dimension aux(num_aux,1-num_ghost:mx+num_ghost)
    double precision :: auxl(num_aux,1-num_ghost:mx+num_ghost), auxr(num_aux,1-num_ghost:mx+num_ghost)
    double precision, intent(out) :: cfl
    integer, intent(in) :: ixy
    integer t
    external rpn2

!f2py intent(in,out) dq1d  
!f2py intent(out) cfl  
!f2py optional dq1d


    if (index_capa.gt.0) then
        dtdx = dt / (dx(ixy)*aux(index_capa,:))
    else
        dtdx = dt/dx(ixy)
    endif
    if (num_dim.gt.1) dq1d=0.d0


    select case(lim_type)
        ! Non-limited reconstruction of components of q (simplest approach)
!        case(0)
!        select case(char_decomp)
!            case(0)
!                call q2qlqr_poly(q1d,ql,qr,mx)
!            case(1)
!                ! wave-based unlimited reconstruction
!                call rpn2(ixy,maxnx,num_eqn,num_waves,num_ghost,mx,&
!                        q1d,q1d,aux,aux,wave,s,amdq,apdq,num_aux)
!                call q2qlqr_poly_wave(q1d,ql,qr,wave,s,mx)
!        end select
        case(1)
        select case(char_decomp)
            case(0)
                ! Fill in TVD reconstruction w/o char. decomp. here
                call tvd2(q1d,ql,qr,mthlim)
            case(1)
                ! wave-based second order reconstruction
                call rpn2(ixy,maxnx,num_eqn,num_waves,num_ghost,mx,&
                        q1d,q1d,aux,aux,wave,s,amdq,apdq,num_aux)
                call tvd2_wave(q1d,ql,qr,wave,s,mthlim)
            case(2)
                ! characteristic-wise second order reconstruction
                call evec(mx,num_eqn,num_ghost,mx,q1d,aux,aux,evl,evr)
                call tvd2_char(q1d,ql,qr,mthlim,evl,evr)
        end select
        case(2)
        select case (char_decomp)
            case (0)
                ! no characteristic decomposition
                call weno_comp(q1d,ql,qr,num_eqn,maxnx,num_ghost)
            case (1)
                ! wave-based reconstruction
                call rpn2(ixy,maxnx,num_eqn,num_waves,num_ghost,mx,&
                        q1d,q1d,aux,aux,wave,s,amdq,apdq,num_aux)
                if (fwave.eqv. .True.) then
                    call weno5_fwave(q1d,ql,qr,wave,s)
                else
                    call weno5_wave(q1d,ql,qr,wave)
                endif
            case (2)
                ! characteristic-wise reconstruction
                call evec(mx,num_eqn,num_ghost,mx,q1d,aux,aux,evl,evr)
                call weno5_char(q1d,ql,qr,evl,evr)
            case (3)
                ! transmission-based reconstruction
                call evec(mx,num_eqn,num_ghost,mx,q1d,aux,aux,evl,evr)
                call weno5_trans(q1d,ql,qr,evl,evr)
            case default
                write(*,*) 'ERROR: Unrecognized characteristic decomposition option'
                write(*,*) 'You should set 0<=char_decomp<=3'
                stop
        end select
        case(3)
            call weno5(q1d,ql,qr,num_eqn,maxnx,num_ghost)
      end select


    ! solve Riemann problem at each interface 
    ! -----------------------------------------
    call rpn2(ixy,maxnx,num_eqn,num_waves,num_ghost,mx,ql,qr,aux,aux, &
              wave,s,amdq,apdq,num_aux)

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
    if (tfluct_solver.eqv. .True.) then
        ! tfluct should be a special solver that uses the parameters aux(i)
        ! to solve a Riemann problem with left state ql(i)
        ! and right state qr(i), and returns a total fluctuation in amdq2
        ! NOTE that here amdq2 is really a total fluctuation (should be
        ! called adq); we do it this way just to avoid declaring more storage
        call tfluct(ixy,maxnx,num_eqn,num_waves,num_ghost,mx,ql,qr, &
                     aux,aux,s,amdq2)

        ! Modify q using fluctuations:
        ! Note this may not correspond to a conservative flux-differencing
        ! for equations not in conservation form.  It is conservative if
        ! adq = f(qr(i)) - f(ql(i)).

        forall (i=1:mx, m=1:num_eqn)
            dq1d(m,i) = dq1d(m,i) - dtdx(i)*(apdq(m,i) + &
                            amdq2(m,i) + amdq(m,i+1))
        end forall

    else
        ! Or we can just swap things around and use the usual Riemann solver
        ! This may be more convenient, but is less efficient. 
        ! For the moment the swapping is done working with the element of the 
        ! vectors qr, ql, auxl, auxr. 
        ! 
        ! TODO: Working with pointers!!!
        
        do i = 1-num_ghost+1,mx+num_ghost
            do m = 1, num_eqn
                qr(m,i-1) = ql(m,i)
                ql(m,i  ) = qr(m,i)
            enddo
        enddo

        if (num_aux .gt. 0) then
             do i = 1-num_ghost+1,mx+num_ghost
                do m = 1, num_aux
                    auxr(m,i-1) = aux(m,i) !aux is not patchdat type
                    auxl(m,i  ) = aux(m,i) !aux is not patchdat type
                enddo
            enddo
        endif
        
        call rpn2(ixy,maxnx,num_eqn,num_waves,num_ghost,mx,ql,qr, &
                 auxl,auxr,wave,s,amdq2,apdq2,num_aux)

        forall(i=1:mx, m=1:num_eqn)
            dq1d(m,i) = dq1d(m,i)-dtdx(i)*(amdq(m,i+1)+ &
                        apdq(m,i)+amdq2(m,i)+apdq2(m,i))
        end forall
    endif

end subroutine flux1
