! ===================================================================
subroutine flux1(q1d,dq1d,ndim,aux,dt,cfl,t,dtdx,ql,qr,wave,s,amdq,apdq, &
                 amdq2,apdq2, &
                 ixy,mcapa,maux,meqn,mx,mwaves,mbc,maxnx,dx,tfluct_solver, &
                 char_decomp,lim_type,mthlim)
! ===================================================================
!
!     # Evaluate (delta t) * dq(t)/dt
!
!     SharpClaw
!     Author: David Ketcheson
!
!     amdq, apdq, amdq2, apdq2, wave, and s are used locally:
!
!     amdq(meqn,1-mbc:mx+mbc) = left-going flux-differences
!     apdq(meqn,1-mbc:mx+mbc) = right-going flux-differences
!        e.g. amdq(m,i) = m'th component of A^- \Delta q from i'th Riemann
!                         problem (between cells i-1 and i).
!
!     wave(1-mbc:mx+mbc, meqn, mwaves) = waves from solution of
!                                           Riemann problems,
!            wave(m,mw,i) = mth component of jump in q across
!                           wave in family mw in Riemann problem between
!                           states i-1 and i.
!
!     s(mwaves,1-mbc:mx+mbc) = wave speeds,
!            s(mw,i) = speed of wave in family mw in Riemann problem between
!                      states i-1 and i.
!
!     t is the time at which we want to evaluate dq/dt, which may not
!      be the current simulation time
! ===================================================================
!
! Add here some info about modifications.
!
! Modified: March 30, 2011
! Authors:  David Ketcheson
!           Matteo Parsani
!
! ===================================================================

    USE reconstruct

    implicit double precision (a-h,o-z)

!f2py intent(in,out) dq1d  
!f2py intent(out) cfl  

    integer :: ndim,mcapa,maux,meqn,mwaves,mbc,maxnx,mx,char_decomp,lim_type
    double precision :: amdq(meqn,1-mbc:mx+mbc),apdq(meqn,1-mbc:mx+mbc)
    double precision :: dtdx(1-mbc:mx+mbc)
    double precision :: amdq2(meqn,1-mbc:mx+mbc),apdq2(meqn,1-mbc:mx+mbc)
    double precision :: qr(meqn,1-mbc:mx+mbc),ql(meqn,1-mbc:mx+mbc)
    double precision :: s(mwaves,1-mbc:mx+mbc),wave(meqn,mwaves,1-mbc:mx+mbc)
    double precision :: q1d(meqn,1-mbc:mx+mbc)
    double precision :: dq1d(meqn,1-mbc:maxnx+mbc)
    dimension aux(maux,1-mbc:mx+mbc)
    double precision :: auxl(maux,1-mbc:mx+mbc), auxr(maux,1-mbc:mx+mbc)
    double precision :: dx(ndim)
    double precision, intent(out) :: cfl
    integer, intent(in) :: ixy, mthlim(mwaves)
    integer t


    if (recon_alloc.eqv..False.) then
        call alloc_recon_workspace(maxnx,mbc,meqn,mwaves,lim_type,char_decomp)
    endif

    if (mcapa.gt.0) then
        dtdx = dt / (dx(ixy)*aux(mcapa,:))
    else
        dtdx = dt/dx(ixy)
    endif
    if (ndim.gt.1) dq1d=0.d0


    select case(lim_type)
        ! Non-limited reconstruction of components of q (simplest approach)
!        case(0)
!        select case(char_decomp)
!            case(0)
!                call q2qlqr_poly(q1d,ql,qr,mx)
!            case(1)
!                ! wave-based unlimited reconstruction
!                call rp1(maxnx,meqn,mwaves,mbc,mx,&
!                        q1d,q1d,aux,aux,wave,s,amdq,apdq)
!                call q2qlqr_poly_wave(q1d,ql,qr,wave,s,mx)
!        end select
!        case(1)
!        select case(char_decomp)
!            case(0)
!                ! Fill in TVD reconstruction w/o char. decomp. here
!                call tvd2(q1d,ql,qr,mthlim)
!            case(1)
!                ! wave-based second order reconstruction
!                call rp1(maxnx,meqn,mwaves,mbc,mx,&
!                        q1d,q1d,aux,aux,wave,s,amdq,apdq)
!                call tvd2_wave(q1d,ql,qr,wave,s,mthlim)
!            case(2)
!                ! characteristic-wise second order reconstruction
!                !call evec(mx,meqn,mbc,mx,q1d,aux,aux,evl,evr)
!                !call tvd2_char(q1d,ql,qr,mthlim,evl,evr)
!        end select
        case(2)
        select case (char_decomp)
            case (0)
                ! no characteristic decomposition
                call weno5(q1d,ql,qr,meqn,maxnx,mbc)
            case (1)
                ! wave-based reconstruction
                !call rp1(maxnx,meqn,mwaves,mbc,mx,&
                !        q1d,q1d,aux,aux,wave,s,amdq,apdq)
                !call weno5_wave(q1d,ql,qr,wave)
                !call weno5_fwave(q1d,ql,qr,wave,s)
            case (2)
                ! characteristic-wise reconstruction
                !call evec(mx,meqn,mbc,mx,q1d,aux,aux,evl,evr)
                !call weno5_char(q1d,ql,qr,evl,evr)
            case (3)
                ! transmission-based reconstruction
                !call evec(mx,meqn,mbc,mx,q1d,aux,aux,evl,evr)
                !call weno5_trans(q1d,ql,qr,evl,evr)
            case default
                write(*,*) 'ERROR: Unrecognized characteristic decomposition option'
                write(*,*) 'You should set 0<=char_decomp<=3'
                stop
        end select
      end select


    ! solve Riemann problem at each interface 
    ! -----------------------------------------
    call rp1(maxnx,meqn,mwaves,mbc,mx,ql,qr,aux,aux, &
              wave,s,amdq,apdq)

    ! compute maximum wave speed:
    cfl = 0.d0
    do mw=1,mwaves
        do i=1,mx+1
            ! if s>0 use dtdx(i) to compute CFL,
            ! if s<0 use dtdx(i-1) to compute CFL:
            cfl = dmax1(cfl, dtdx(i)*s(mw,i), -dtdx(i-1)*s(mw,i))
        enddo
    enddo

    ! Find total fluctuation within each cell
!    if (tfluct_solver==1) then
!        ! tfluct should be a special solver that uses the parameters aux(i)
!        ! to solve a Riemann problem with left state ql(i)
!        ! and right state qr(i), and returns a total fluctuation in amdq2
!        ! NOTE that here amdq2 is really a total fluctuation (should be
!        ! called adq); we do it this way just to avoid declaring more storage
!        call tfluct(ixy,maxnx,meqn,mwaves,mbc,mx,ql,qr, &
!                     aux,aux,s,amdq2)
!
!        ! Modify q using fluctuations:
!        ! Note this may not correspond to a conservative flux-differencing
!        ! for equations not in conservation form.  It is conservative if
!        ! adq = f(qr(i)) - f(ql(i)).
!
!        forall (i=1:mx, m=1:meqn)
!            dq1d(m,i) = dq1d(m,i) - dtdx(i)*(apdq(m,i) + &
!                            amdq2(m,i) + amdq(m,i+1))
!        end forall
!
!    else
        ! Or we can just swap things around and use the usual Riemann solver
        ! This may be more convenient, but is less efficient. 
        ! For the moment the swapping is done working with the element of the 
        ! vectors qr, ql, auxl, auxr. 
        ! 
        ! TODO: Working with pointers!!!
        
        do i = 1-mbc+1,mx+mbc
            do m = 1, meqn
                qr(m,i-1) = ql(m,i)
                ql(m,i  ) = qr(m,i)
            enddo
        enddo

        if (maux .gt. 0) then
             do i = 1-mbc+1,mx+mbc
                do m = 1, maux
                    auxr(m,i-1) = aux(m,i) !aux is not griddat type
                    auxl(m,i  ) = aux(m,i) !aux is not griddat type
                enddo
            enddo
        endif
        
        call rp1(maxnx,meqn,mwaves,mbc,mx,ql,qr, &
                 auxl,auxr,wave,s,amdq2,apdq2)

        forall(i=1:mx, m=1:meqn)
            dq1d(m,i) = dq1d(m,i)-dtdx(i)*(amdq(m,i+1)+ &
                        apdq(m,i)+amdq2(m,i)+apdq2(m,i))
        end forall
!    endif

end subroutine flux1
