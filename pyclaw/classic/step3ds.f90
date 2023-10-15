! ==================================================================
subroutine step3ds(maxm,num_eqn,num_waves,num_ghost,mx,my, &
                   mz,qold,qnew,aux,dx,dy,dz,dt,method,mthlim,cfl, &
                   qadd,fadd,gadd,hadd,q1d,dtdx1d,dtdy1d,dtdz1d, &
                   aux1,aux2,aux3,num_aux,work,mwork,idir,use_fwave,rpn3,rpt3,rptt3)
! ==================================================================

!  Take one time step, updating q, to be used with
!  dimensional splitting.
!  On entry, qold and qnew should be identical and give the
!     initial data for this step.
!  On exit, qnew returns values at the end of the time step.
!     qold is unchanged.

!  This is a simplified version of the subroutine step3d
!  Since only qadd and fadd is used in the flux updates,
!  the solution updates below are considerably simplified
!  compared to step3d.

! -----------------------------------------------------------------
!  NOTE! Since dimensional splitting is used, it is possible
!        to reduce the memory requirement, i.e. the size of
!        the work array. It could be reduced with
! 
!        (maxm + 2*num_ghost)*(37*num_eqn + 6*num_aux),
! 
!        when also possible reductions in flux3 are included.
!        However, this term is small compared to the dominating
!        term (mx+2num_ghost)(my+2mb)*(mz+2num_ghost).
! -----------------------------------------------------------------

    implicit real*8(a-h,o-z)
    external rpn3,rpt3,rptt3
    dimension qold(num_eqn, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
    dimension qnew(num_eqn, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
    dimension  q1d(num_eqn,1-num_ghost:maxm+num_ghost)
    dimension qadd(num_eqn,1-num_ghost:maxm+num_ghost)
    dimension fadd(num_eqn,1-num_ghost:maxm+num_ghost)
    dimension gadd(num_eqn,2,-1:1,1-num_ghost:maxm+num_ghost)
    dimension hadd(num_eqn,2,-1:1,1-num_ghost:maxm+num_ghost)
    dimension aux(num_aux, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
    dimension aux1(num_aux,1-num_ghost:maxm+num_ghost,3)
    dimension aux2(num_aux,1-num_ghost:maxm+num_ghost,3)
    dimension aux3(num_aux,1-num_ghost:maxm+num_ghost,3)
    dimension dtdx1d(1-num_ghost:maxm+num_ghost)
    dimension dtdy1d(1-num_ghost:maxm+num_ghost)
    dimension dtdz1d(1-num_ghost:maxm+num_ghost)
    dimension method(7),mthlim(num_waves)
    dimension work(mwork)
    logical ::          use_fwave

    common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom

!f2py intent(out) cfl
!f2py intent(in,out) qnew
!f2py optional q1d, qadd, fadd, gadd, hadd, dtdx1d, dtdy1d, dtdz1d

! Dummy interfaces just so f2py doesn't complain:
!f2py real(DP) x
!f2py x=rpn3(x)
!f2py x=rpt3(x)
!f2py x=rptt3(x)

    ! partition work array into pieces needed for local storage in
    ! flux2 routine.  Find starting index of each piece:

    i0wave     = 1
    i0s        = i0wave     + (maxm+2*num_ghost)*num_eqn*num_waves
    i0amdq     = i0s        + (maxm+2*num_ghost)*num_waves
    i0apdq     = i0amdq     + (maxm+2*num_ghost)*num_eqn
    i0cqxx     = i0apdq     + (maxm+2*num_ghost)*num_eqn
    i0bmamdq   = i0cqxx     + (maxm+2*num_ghost)*num_eqn
    i0bmapdq   = i0bmamdq   + (maxm+2*num_ghost)*num_eqn
    i0bpamdq   = i0bmapdq   + (maxm+2*num_ghost)*num_eqn
    i0bpapdq   = i0bpamdq   + (maxm+2*num_ghost)*num_eqn
    i0cmamdq   = i0bpapdq   + (maxm+2*num_ghost)*num_eqn
    i0cmapdq   = i0cmamdq   + (maxm+2*num_ghost)*num_eqn
    i0cpamdq   = i0cmapdq   + (maxm+2*num_ghost)*num_eqn
    i0cpapdq   = i0cpamdq   + (maxm+2*num_ghost)*num_eqn
    i0cmamdq2  = i0cpapdq   + (maxm+2*num_ghost)*num_eqn
    i0cmapdq2  = i0cmamdq2  + (maxm+2*num_ghost)*num_eqn
    i0cpamdq2  = i0cmapdq2  + (maxm+2*num_ghost)*num_eqn
    i0cpapdq2  = i0cpamdq2  + (maxm+2*num_ghost)*num_eqn
    i0bmcqxxp  = i0cpapdq2  + (maxm+2*num_ghost)*num_eqn
    i0bmcqxxm  = i0bmcqxxp  + (maxm+2*num_ghost)*num_eqn
    i0bpcqxxp  = i0bmcqxxm   + (maxm+2*num_ghost)*num_eqn
    i0bpcqxxm  = i0bpcqxxp   + (maxm+2*num_ghost)*num_eqn
    i0cmcqxxp  = i0bpcqxxm   + (maxm+2*num_ghost)*num_eqn
    i0cmcqxxm  = i0cmcqxxp   + (maxm+2*num_ghost)*num_eqn
    i0cpcqxxp  = i0cmcqxxm   + (maxm+2*num_ghost)*num_eqn
    i0cpcqxxm  = i0cpcqxxp   + (maxm+2*num_ghost)*num_eqn
    i0bmcmamdq = i0cpcqxxm   + (maxm+2*num_ghost)*num_eqn
    i0bmcmapdq = i0bmcmamdq + (maxm+2*num_ghost)*num_eqn
    i0bpcmamdq = i0bmcmapdq + (maxm+2*num_ghost)*num_eqn
    i0bpcmapdq = i0bpcmamdq + (maxm+2*num_ghost)*num_eqn
    i0bmcpamdq = i0bpcmapdq + (maxm+2*num_ghost)*num_eqn
    i0bmcpapdq = i0bmcpamdq + (maxm+2*num_ghost)*num_eqn
    i0bpcpamdq = i0bmcpapdq + (maxm+2*num_ghost)*num_eqn
    i0bpcpapdq = i0bpcpamdq + (maxm+2*num_ghost)*num_eqn
    iused      = i0bpcpapdq + (maxm+2*num_ghost)*num_eqn - 1

    if (iused > mwork) then
        ! This shouldn't happen due to checks in claw2
        write(6,*) '*** not enough work space in step2'
        write(6,*) '*** iused = ', iused, '   mwork =',mwork
        stop
    endif

    index_capa = method(6)
    ! num_aux = method(7)
    cfl = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz

    if (index_capa == 0) then
        ! no capa array:
        do i=1-num_ghost,maxm+num_ghost
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
            dtdz1d(i) = dtdz
        end do
    endif


    if( idir == 1)then
    
        ! perform x-sweeps
        ! ==================
    
        do k = 0,mz+1
            do j = 0,my+1
            
                forall (m = 1:num_eqn, i = 1-num_ghost:mx+num_ghost)
                    ! copy data along a slice into 1d array:
                    q1d(m,i) = qold(m,i,j,k)
                end forall
            
                if (index_capa > 0)  then
                    do i = 1-num_ghost, mx+num_ghost
                        dtdx1d(i) = dtdx / aux(index_capa,i,j,k)
                    end do
                endif
            
                ! Since dimensional splitting is used, only aux2 is needed.
            
                if (num_aux > 0)  then
                    forall (ma=1:num_aux,i= 1-num_ghost:mx+num_ghost, ka = -1:1)
                    aux2(ma,i,2+ka) = aux(ma,i,j,k+ka)
                    end forall
                endif
            
                ! Store the value of j and k along this slice in the common block
                ! comxyt in case it is needed in the Riemann solver (for
                ! variable coefficient problems)
            
                jcom = j
                kcom = k
            
                ! compute modifications qadd and fadd along this slice
            
                call flux3(1,maxm,num_eqn,num_waves,num_ghost,mx, &
                           q1d,dtdx1d,dtdy,dtdz,dummy1,aux2,dummy3,num_aux, &
                           method,mthlim,qadd,fadd,gadd,hadd,cfl1d, &
                           work(i0wave),work(i0s),work(i0amdq), &
                           work(i0apdq),work(i0cqxx), &
                           work(i0bmamdq),work(i0bmapdq), &
                           work(i0bpamdq),work(i0bpapdq), &
                           work(i0cmamdq),work(i0cmapdq), &
                           work(i0cpamdq),work(i0cpapdq), &
                           work(i0cmamdq2),work(i0cmapdq2), &
                           work(i0cpamdq2),work(i0cpapdq2), &
                           work(i0bmcqxxp),work(i0bpcqxxp), &
                           work(i0bmcqxxm),work(i0bpcqxxm), &
                           work(i0cmcqxxp),work(i0cpcqxxp), &
                           work(i0cmcqxxm),work(i0cpcqxxm), &
                           work(i0bmcmamdq),work(i0bmcmapdq), &
                           work(i0bpcmamdq),work(i0bpcmapdq), &
                           work(i0bmcpamdq),work(i0bmcpapdq), &
                           work(i0bpcpamdq),work(i0bpcpapdq), &
                           use_fwave,rpn3,rpt3,rptt3)
                        
                cfl = dmax1(cfl,cfl1d)
            
                ! update qnew by flux differencing.
                ! (rather than maintaining arrays f, g and h for the total fluxes,
                ! the modifications are used immediately to update qnew
                ! in order to save storage.)
            
                if(index_capa == 0)then
                    ! no capa array.  Standard flux differencing:
                    forall (m = 1:num_eqn, i = 1:mx)
                        qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,i) &
                            - dtdx * (fadd(m,i+1) - fadd(m,i))
                    end forall
                else
                    ! with capa array
                    forall (m = 1:num_eqn, i = 1:mx)
                        qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,i) &
                            - dtdx * (fadd(m,i+1) - fadd(m,i)) &
                            / aux(index_capa,i,j,k)
                    end forall
                endif
            end do
        end do
    
    else if( idir == 2 )then
    
        ! perform y sweeps
        ! ==================
    
        do k = 0, mz+1
            do i = 0, mx+1
            
                forall (m = 1:num_eqn, j = 1-num_ghost:my+num_ghost)
                    ! copy data along a slice into 1d array:
                    q1d(m,j) = qold(m,i,j,k)
                end forall
            
                if (index_capa > 0)  then
                    do j = 1-num_ghost, my+num_ghost
                        dtdy1d(j) = dtdy / aux(index_capa,i,j,k)
                    end do
                endif
            
                ! Since dimensional splitting is used, only aux2 is needed.
            
                if (num_aux > 0)  then
                    ! There's a decent chance aux2 will all fit into cache,
                    ! so keep accesses to aux as contiguous as possible.
                    forall (ma=1:num_aux,ia= -1:1, j = 1-num_ghost:my+num_ghost)
                    aux2(ma, j, 2+ia) = aux(ma, i+ia, j, k)
                    end forall
                endif
            
                ! Store the value of i and k along this slice in the common block
                ! comxyzt in case it is needed in the Riemann solver (for
                ! variable coefficient problems)
            
                icom = i
                kcom = k
            
                ! compute modifications qadd and fadd along this slice
            
                call flux3(2,maxm,num_eqn,num_waves,num_ghost,my, &
                           q1d,dtdy1d,dtdz,dtdx,dummy1,aux2,dummy3,num_aux, &
                           method,mthlim,qadd,fadd,gadd,hadd,cfl1d, &
                           work(i0wave),work(i0s),work(i0amdq), &
                           work(i0apdq),work(i0cqxx), &
                           work(i0bmamdq),work(i0bmapdq), &
                           work(i0bpamdq),work(i0bpapdq), &
                           work(i0cmamdq),work(i0cmapdq), &
                           work(i0cpamdq),work(i0cpapdq), &
                           work(i0cmamdq2),work(i0cmapdq2), &
                           work(i0cpamdq2),work(i0cpapdq2), &
                           work(i0bmcqxxp),work(i0bpcqxxp), &
                           work(i0bmcqxxm),work(i0bpcqxxm), &
                           work(i0cmcqxxp),work(i0cpcqxxp), &
                           work(i0cmcqxxm),work(i0cpcqxxm), &
                           work(i0bmcmamdq),work(i0bmcmapdq), &
                           work(i0bpcmamdq),work(i0bpcmapdq), &
                           work(i0bmcpamdq),work(i0bmcpapdq), &
                           work(i0bpcpamdq),work(i0bpcpapdq), &
                           use_fwave,rpn3,rpt3,rptt3)
            
                cfl = dmax1(cfl,cfl1d)
            
                ! update qnew by flux differencing.
                ! Note that the roles of the flux updates are changed.
                ! fadd - modifies the g-fluxes
            
                if( index_capa == 0)then
                    ! no capa array.  Standard flux differencing:
                    forall (m = 1:num_eqn, j = 1:my)
                        qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,j) &
                            - dtdy * (fadd(m,j+1) - fadd(m,j))
                    end forall
                else
                    ! with capa array.
                    forall (m = 1:num_eqn, j = 1:my)
                        qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,j) &
                            - dtdy * (fadd(m,j+1) - fadd(m,j)) &
                            / aux(index_capa,i,j,k)
                    end forall
                endif
            end do
        end do
    
    else
    
        ! perform z sweeps
        ! ==================
    
    
        do j = 0, my+1
            do i = 0, mx+1
            
                forall (m = 1:num_eqn, k = 1-num_ghost:mz+num_ghost)
                    ! copy data along a slice into 1d array:
                    q1d(m,k) = qold(m,i,j,k)
                end forall
            
                if (index_capa > 0)  then
                    do k = 1-num_ghost, mz+num_ghost
                        dtdz1d(k) = dtdz / aux(index_capa,i,j,k)
                    end do
                endif
            
                ! Since dimensional splitting is used, only aux2 is needed.
            
                if (num_aux > 0)  then
                    ! There's a decent chance aux2 will all fit into cache,
                    ! so keep accesses to aux as contiguous as possible.
                    forall (ma=1:num_aux,ja= -1:1, k = 1-num_ghost:mz+num_ghost)
                        aux2(ma, k, 2+ja) = aux(ma, i, j+ja, k)
                    end forall
                endif
            
                ! Store the value of i and j along this slice in the common block
                ! comxyzt in case it is needed in the Riemann solver (for
                ! variable coefficient problems)
            
                icom = i
                jcom = j
            
                ! compute modifications qadd and fadd along this slice
            
                call flux3(3,maxm,num_eqn,num_waves,num_ghost,mz, &
                           q1d,dtdz1d,dtdx,dtdy,dummy1,aux2,dummy3,num_aux, &
                           method,mthlim,qadd,fadd,gadd,hadd,cfl1d, &
                           work(i0wave),work(i0s),work(i0amdq), &
                           work(i0apdq),work(i0cqxx), &
                           work(i0bmamdq),work(i0bmapdq), &
                           work(i0bpamdq),work(i0bpapdq), &
                           work(i0cmamdq),work(i0cmapdq), &
                           work(i0cpamdq),work(i0cpapdq), &
                           work(i0cmamdq2),work(i0cmapdq2), &
                           work(i0cpamdq2),work(i0cpapdq2), &
                           work(i0bmcqxxp),work(i0bpcqxxp), &
                           work(i0bmcqxxm),work(i0bpcqxxm), &
                           work(i0cmcqxxp),work(i0cpcqxxp), &
                           work(i0cmcqxxm),work(i0cpcqxxm), &
                           work(i0bmcmamdq),work(i0bmcmapdq), &
                           work(i0bpcmamdq),work(i0bpcmapdq), &
                           work(i0bmcpamdq),work(i0bmcpapdq), &
                           work(i0bpcpamdq),work(i0bpcpapdq), &
                           use_fwave,rpn3,rpt3,rptt3)
            
                cfl = dmax1(cfl,cfl1d)
            
                ! update qnew by flux differencing.
                ! Note that the roles of the flux updates are changed.
                ! fadd - modifies the h-fluxes
            
                if(index_capa == 0)then
                    ! no capa array. Standard flux differencing:
                    forall (m = 1:num_eqn, k = 1:mz)
                        qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,k) &
                            - dtdz * (fadd(m,k+1) - fadd(m,k))
                    end forall
                else
                    ! with capa array
                    forall (m = 1:num_eqn, k = 1:mz)
                        qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,k) &
                            - dtdz * (fadd(m,k+1) - fadd(m,k)) &
                            / aux(index_capa,i,j,k)
                    end forall
                endif
            end do
        end do
    
    endif

    return
end subroutine step3ds


