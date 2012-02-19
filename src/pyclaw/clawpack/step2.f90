!     ==========================================================
    subroutine step2(maxm,num_eqn,num_waves,num_aux,num_ghost,mx,my, &
    qold,qnew,aux,dx,dy,dt,method,mthlim,cfl, &
    qadd,fadd,gadd,q1d,dtdx1d,dtdy1d, &
    aux1,aux2,aux3,work,mwork,use_fwave,rpn2,rpt2)
!     ==========================================================

!     # Take one time step, updating q.
!     # On entry, qold and qnew should be identical and give the
!     #    initial data for this step
!     # On exit, qnew returns values at the end of the time step.
!     #    qold is unchanged.

!     # qadd is used to return increments to q from flux2
!     # fadd and gadd are used to return flux increments from flux2.
!     # See the flux2 documentation for more information.


    implicit double precision (a-h,o-z)
    external :: rpn2,rpt2
    double precision :: qold(num_eqn,1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision :: qnew(num_eqn,1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision ::  q1d(num_eqn,1-num_ghost:maxm+num_ghost)
    double precision :: qadd(num_eqn,1-num_ghost:maxm+num_ghost)
    double precision :: fadd(num_eqn,1-num_ghost:maxm+num_ghost)
    double precision :: gadd(num_eqn, 2, 1-num_ghost:maxm+num_ghost)
    double precision :: aux(num_aux, 1-num_ghost:mx+num_ghost, &
    1-num_ghost:my+num_ghost)
    double precision :: aux1(num_aux, 1-num_ghost:maxm+num_ghost)
    double precision :: aux2(num_aux, 1-num_ghost:maxm+num_ghost)
    double precision :: aux3(num_aux, 1-num_ghost:maxm+num_ghost)

    double precision :: dtdx1d(1-num_ghost:maxm+num_ghost)
    double precision :: dtdy1d(1-num_ghost:maxm+num_ghost)
    integer :: method(7),mthlim(num_waves)
    logical ::          use_fwave
    double precision :: work(mwork)

    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

!f2py intent(out) cfl
!f2py intent(in,out) qnew
!f2py optional q1d, qadd, fadd, gadd, dtdx1d, dtdy1d

! Dummy interfaces just so f2py doesn't complain:
!f2py real(DP) x
!f2py x=rpn2(x)
!f2py x=rpt2(x)



!     # partition work array into pieces needed for local storage in
!     # flux2 routine.  Find starting index of each piece:

    i0wave = 1
    i0s = i0wave + (maxm+2*num_ghost)*num_eqn*num_waves
    i0amdq = i0s + (maxm+2*num_ghost)*num_waves
    i0apdq = i0amdq + (maxm+2*num_ghost)*num_eqn
    i0cqxx = i0apdq + (maxm+2*num_ghost)*num_eqn
    i0bmadq = i0cqxx + (maxm+2*num_ghost)*num_eqn
    i0bpadq = i0bmadq + (maxm+2*num_ghost)*num_eqn
    iused = i0bpadq + (maxm+2*num_ghost)*num_eqn - 1

    if (iused > mwork) then
    !        # This shouldn't happen due to checks in claw2
        write(6,*) '*** not enough work space in step2'
        write(6,*) '*** iused = ', iused, '   mwork =',mwork
        stop
    endif


    index_capa = method(6)
    num_aux = method(7)
    cfl = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy

    if (index_capa == 0) then
    !        # no capa array:
        do 5 i=1-num_ghost,maxm+num_ghost
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
        5 END DO
    endif


!     # perform x-sweeps
!     ==================

    do 50 j = 0,my+1
    
    !        # copy data along a slice into 1d arrays:
        do 21 m=1,num_eqn
            do 20 i = 1-num_ghost, mx+num_ghost
                q1d(m,i) = qold(m,i,j)
            20 END DO
        21 END DO
    
        if (index_capa > 0)  then
            do 22 i = 1-num_ghost, mx+num_ghost
                dtdx1d(i) = dtdx/aux(index_capa,i,j)
            22 END DO
        endif
    
        if (num_aux > 0)  then
            do 25 ma=1,num_aux
                do 24 i = 1-num_ghost, mx+num_ghost
                    aux1(ma,i) = aux(ma,i,j-1)
                    aux2(ma,i) = aux(ma,i,j  )
                    aux3(ma,i) = aux(ma,i,j+1)
                24 END DO
            25 END DO
        endif
    
    !     # Store the value of j along this slice in the common block
    !        # comxyt in case it is needed in the Riemann solver (for
    !        # variable coefficient problems)
        jcom = j
    
    !        # compute modifications fadd and gadd to fluxes along this slice:
        call flux2(1,maxm,num_eqn,num_waves,num_aux,num_ghost,mx, &
        q1d,dtdx1d,aux1,aux2,aux3,method,mthlim, &
        qadd,fadd,gadd,cfl1d, &
        work(i0wave),work(i0s),work(i0amdq),work(i0apdq), &
        work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2, &
        use_fwave)
        cfl = dmax1(cfl,cfl1d)
    
    !        # update qnew by flux differencing.
    !        # (rather than maintaining arrays f and g for the total fluxes,
    !        # the modifications are used immediately to update qnew
    !        # in order to save storage.)
    
        if (index_capa == 0) then
        
        !            # no capa array.  Standard flux differencing:
            do 31 m=1,num_eqn
                do 30 i=1,mx
                    qnew(m,i,j) = qnew(m,i,j) + qadd(m,i) &
                    - dtdx * (fadd(m,i+1) - fadd(m,i)) &
                    - dtdy * (gadd(m,2,i) - gadd(m,1,i))
                    qnew(m,i,j-1) = qnew(m,i,j-1) - dtdy * gadd(m,1,i)
                    qnew(m,i,j+1) = qnew(m,i,j+1) + dtdy * gadd(m,2,i)
                30 END DO
            31 END DO
        
        else
        
        !            # with capa array.
            do 41 m=1,num_eqn
                do 40 i=1,mx
                    qnew(m,i,j) = qnew(m,i,j) + qadd(m,i) &
                    - (dtdx * (fadd(m,i+1) - fadd(m,i)) &
                    +  dtdy * (gadd(m,2,i) - gadd(m,1,i))) &
                    / aux(index_capa,i,j)
                    qnew(m,i,j-1) = qnew(m,i,j-1) - dtdy * gadd(m,1,i) &
                    / aux(index_capa,i,j-1)
                    qnew(m,i,j+1) = qnew(m,i,j+1) + dtdy * gadd(m,2,i) &
                    / aux(index_capa,i,j+1)
                40 END DO
            41 END DO
        endif
    50 END DO



!     # perform y sweeps
!     ==================


    do 100 i = 0, mx+1
    
    !        # copy data along a slice into 1d arrays:
        do 71 m=1,num_eqn
            do 70 j = 1-num_ghost, my+num_ghost
                q1d(m,j) = qold(m,i,j)
            70 END DO
        71 END DO
    
        if (index_capa > 0)  then
            do 72 j = 1-num_ghost, my+num_ghost
                dtdy1d(j) = dtdy / aux(index_capa,i,j)
            72 END DO
        endif
    
        if (num_aux > 0)  then
            do 75 ma=1,num_aux
                do 74 j = 1-num_ghost, my+num_ghost
                    aux1(ma,j) = aux(ma,i-1,j)
                    aux2(ma,j) = aux(ma,i,  j)
                    aux3(ma,j) = aux(ma,i+1,j)
                74 END DO
            75 END DO
        endif
    
    
    !     # Store the value of i along this slice in the common block
    !        # comxyt in case it is needed in the Riemann solver (for
    !        # variable coefficient problems)
        icom = i
    
    !        # compute modifications fadd and gadd to fluxes along this slice:
        call flux2(2,maxm,num_eqn,num_waves,num_aux,num_ghost,my, &
        q1d,dtdy1d,aux1,aux2,aux3,method,mthlim, &
        qadd,fadd,gadd,cfl1d, &
        work(i0wave),work(i0s),work(i0amdq),work(i0apdq), &
        work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2, &
        use_fwave)
    
        cfl = dmax1(cfl,cfl1d)
    
    !        # update qnew by flux differencing.
    !        # Note that the roles of fadd and gadd are reversed for
    !        # the y-sweeps -- fadd is the modification to g-fluxes and
    !        # gadd is the modification to f-fluxes to the left and right.
    
        if (index_capa == 0) then
        
        !            # no capa array.  Standard flux differencing:
            do 81 m=1,num_eqn
                do 80 j=1,my
                    qnew(m,i,j) = qnew(m,i,j) + (qadd(m,j) &
                    - dtdy * (fadd(m,j+1) - fadd(m,j)) &
                    - dtdx * (gadd(m,2,j) - gadd(m,1,j)))
                    qnew(m,i-1,j) = qnew(m,i-1,j) - dtdx * gadd(m,1,j)
                    qnew(m,i+1,j) = qnew(m,i+1,j) + dtdx * gadd(m,2,j)
                80 END DO
            81 END DO
        
        else
        
        !            # with capa array.
            do 91 m=1,num_eqn
                do 90 j=1,my
                    qnew(m,i,j) = qnew(m,i,j) + qadd(m,j) &
                    - (dtdy * (fadd(m,j+1) - fadd(m,j)) &
                    + dtdx * (gadd(m,2,j) - gadd(m,1,j))) &
                    / aux(index_capa,i,j)
                    qnew(m,i-1,j) = qnew(m,i-1,j) - dtdx * gadd(m,1,j) &
                    / aux(index_capa,i-1,j)
                    qnew(m,i+1,j) = qnew(m,i+1,j) + dtdx * gadd(m,2,j) &
                    / aux(index_capa,i+1,j)
                90 END DO
            91 END DO
        endif
    100 END DO

    return
    end subroutine step2
