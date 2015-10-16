
! ==================================================================
subroutine step3ds(maxm,nthreads,num_eqn,num_waves,num_ghost,mx,my, &
     mz,qold,qnew,aux,dx,dy,dz,dt,method,mthlim,cfl, &
     qadd,fadd,gadd,hadd,q1d,dtdx1d,dtdy1d,dtdz1d,aux1,aux2,aux3,&
     num_aux,work,mwork,idir,use_fwave,rpn3,rpt3,rptt3)
! ==================================================================

  ! Take one time step, updating q, to be used with
  ! dimensional splitting.
  ! On entry, qold and qnew should be identical and give the
  !    initial data for this step.
  ! On exit, qnew returns values at the end of the time step.
  !    qold is unchanged.
  !
  ! This is a simplified version of the subroutine step3d
  ! Since only qadd and fadd is used in the flux updates,
  ! the solution updates below are considerably simplified
  ! compared to step3d.
  
  !-----------------------------------------------------------------
  ! NOTE! Since dimensional splitting is used, it is possible
  !       to reduce the memory requirement, i.e. the size of
  !       the work array. It could be reduced with
  !
  !       (maxm + 2*num_ghost)*(37*num_eqn + 6*num_aux),
  !
  !       when also possible reductions in flux3 are included.
  !       However, this term is small compared to the dominating
  !       term (mx+2num_ghost)(my+2mb)*(mz+2num_ghost).
  !-----------------------------------------------------------------
  
  !$ use omp_lib

  implicit real*8(a-h,o-z)
  external rpn3,rpt3,rptt3
  dimension qold(num_eqn, 1-num_ghost:mx+num_ghost, &
       1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
  dimension qnew(num_eqn, 1-num_ghost:mx+num_ghost, &
       1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
  dimension  q1d(num_eqn,1-num_ghost:maxm+num_ghost,mx+nthreads-mx) ! Adding and subtracting mx is a hack to trick f2py, so it does not make nthreads optional.
  dimension qadd(num_eqn,1-num_ghost:maxm+num_ghost,mx+nthreads-mx)
  dimension fadd(num_eqn,1-num_ghost:maxm+num_ghost,mx+nthreads-mx)
  dimension gadd(num_eqn,2,-1:1,1-num_ghost:maxm+num_ghost,mx+nthreads-mx)
  dimension hadd(num_eqn,2,-1:1,1-num_ghost:maxm+num_ghost,mx+nthreads-mx)
  dimension aux(num_aux, 1-num_ghost:mx+num_ghost, &
       1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
  dimension aux1(num_aux,1-num_ghost:maxm+num_ghost,3,mx+nthreads-mx)
  dimension aux2(num_aux,1-num_ghost:maxm+num_ghost,3,mx+nthreads-mx)
  dimension aux3(num_aux,1-num_ghost:maxm+num_ghost,3,mx+nthreads-mx)
  dimension dtdx1d(1-num_ghost:maxm+num_ghost,mx+nthreads-mx)
  dimension dtdy1d(1-num_ghost:maxm+num_ghost,mx+nthreads-mx)
  dimension dtdz1d(1-num_ghost:maxm+num_ghost,mx+nthreads-mx)
  dimension method(7),mthlim(num_waves)
  dimension work(mwork)
  logical :: use_fwave

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
  ! flux3 routine.  Find starting index of each piece:

  nthreads=1 ! Serial
  !$OMP PARALLEL
  !$OMP SINGLE
  !$ nthreads = omp_get_num_threads()
  !$OMP END SINGLE
  !$OMP END PARALLEL

  nsiz = (maxm+2*num_ghost)*num_eqn
  nsiz_w = nsiz * num_waves
  nsiz_s = (maxm+2*num_ghost)*num_waves

  i0wave     = 1
  i0s        = i0wave     + nsiz_w*nthreads
  i0amdq     = i0s        + nsiz_s*nthreads
  i0apdq     = i0amdq     + nsiz*nthreads
  i0cqxx     = i0apdq     + nsiz*nthreads
  i0bmamdq   = i0cqxx     + nsiz*nthreads
  i0bmapdq   = i0bmamdq   + nsiz*nthreads
  i0bpamdq   = i0bmapdq   + nsiz*nthreads
  i0bpapdq   = i0bpamdq   + nsiz*nthreads
  i0cmamdq   = i0bpapdq   + nsiz*nthreads
  i0cmapdq   = i0cmamdq   + nsiz*nthreads
  i0cpamdq   = i0cmapdq   + nsiz*nthreads
  i0cpapdq   = i0cpamdq   + nsiz*nthreads
  i0cmamdq2  = i0cpapdq   + nsiz*nthreads
  i0cmapdq2  = i0cmamdq2  + nsiz*nthreads
  i0cpamdq2  = i0cmapdq2  + nsiz*nthreads
  i0cpapdq2  = i0cpamdq2  + nsiz*nthreads
  i0bmcqxxp  = i0cpapdq2  + nsiz*nthreads
  i0bmcqxxm  = i0bmcqxxp  + nsiz*nthreads
  i0bpcqxxp  = i0bmcqxxm  + nsiz*nthreads
  i0bpcqxxm  = i0bpcqxxp  + nsiz*nthreads
  i0cmcqxxp  = i0bpcqxxm  + nsiz*nthreads
  i0cmcqxxm  = i0cmcqxxp  + nsiz*nthreads
  i0cpcqxxp  = i0cmcqxxm  + nsiz*nthreads
  i0cpcqxxm  = i0cpcqxxp  + nsiz*nthreads
  i0bmcmamdq = i0cpcqxxm  + nsiz*nthreads
  i0bmcmapdq = i0bmcmamdq + nsiz*nthreads
  i0bpcmamdq = i0bmcmapdq + nsiz*nthreads
  i0bpcmapdq = i0bpcmamdq + nsiz*nthreads
  i0bmcpamdq = i0bpcmapdq + nsiz*nthreads
  i0bmcpapdq = i0bmcpamdq + nsiz*nthreads
  i0bpcpamdq = i0bmcpapdq + nsiz*nthreads
  i0bpcpapdq = i0bpcpamdq + nsiz*nthreads
  iused      = i0bpcpapdq + nsiz*nthreads - 1

  if (iused > mwork) then
     ! This shouldn't happen due to checks in claw3
     write(6,*) '*** not enough work space in step3'
     write(6,*) '*** iused = ', iused, '   mwork =',mwork
     stop
  endif

  index_capa = method(6)
  !num_aux = method(7)
  cfl = 0.d0
  dtdx = dt/dx
  dtdy = dt/dy
  dtdz = dt/dz
  
  if (index_capa == 0) then
     ! no capa array:
     !$OMP PARALLEL PRIVATE(i,me,me1)
     me = 0
     !$ me = omp_get_thread_num()
     me1 = me + 1
     do i=1-num_ghost,maxm+num_ghost
        dtdx1d(i,me1) = dtdx
        dtdy1d(i,me1) = dtdy
        dtdz1d(i,me1) = dtdz
     end do
     !$OMP END PARALLEL
  endif

  if( idir == 1)then
    
     ! perform x-sweeps
     ! ==================
     !$OMP PARALLEL PRIVATE(me,me1,m,i,ma,ka,cfl1d) REDUCTION(MAX:cfl)
     me = 0
     !$ me = omp_get_thread_num()
     me1 = me + 1
     cfl1d = 0.d0

     ! Guided schedule seems to produce a small performance increase
     ! relative to static or dynamic for problems with uniform work
     ! per column.  For problems with very nonuniform work per
     ! column, either guided or dynamic is necessary to keep all the
     ! threads busy.
     
     !$OMP DO COLLAPSE(2) SCHEDULE(GUIDED)
     do k = 0,mz+1
        do j = 0,my+1
           
           forall (m = 1:num_eqn, i = 1-num_ghost:mx+num_ghost)
              ! copy data along a slice into 1d array:
              q1d(m,i,me1) = qold(m,i,j,k)
           end forall
           
           if (index_capa > 0) then
              do i = 1-num_ghost, mx+num_ghost
                 dtdx1d(i,me1) = dtdx / aux(index_capa,i,j,k)
              end do
           endif
            
           ! Since dimensional splitting is used, only aux2 is needed.
           if (num_aux > 0)  then
              forall (ma=1:num_aux,i= 1-num_ghost:mx+num_ghost, ka = -1:1)
                 aux2(ma,i,2+ka,me1) = aux(ma,i,j,k+ka)
              end forall
           endif
            
           ! Store the value of j and k along this slice in the common block
           ! comxyt in case it is needed in the Riemann solver (for
           ! variable coefficient problems)
           
           jcom = j
           kcom = k
            
           ! compute modifications qadd and fadd along this slice
            
           call flux3(1,maxm,num_eqn,num_waves,num_ghost,mx, &
                q1d(1,1-num_ghost,me1),dtdx1d(1-num_ghost,me1),dtdy,dtdz,dummy1,&
                aux2(1,1-num_ghost,1,me1),dummy3,num_aux, &
                method,mthlim,qadd(1,1-num_ghost,me1),fadd(1,1-num_ghost,me1),&
                gadd(1,1,-1,1-num_ghost,me1),hadd(1,1,-1,1-num_ghost,me1),cfl1d, &
                work(i0wave+me*nsiz_w),work(i0s+me*nsiz_s),work(i0amdq+me*nsiz), &
                work(i0apdq+me*nsiz),work(i0cqxx+me*nsiz), &
                work(i0bmamdq+me*nsiz),work(i0bmapdq+me*nsiz), &
                work(i0bpamdq+me*nsiz),work(i0bpapdq+me*nsiz), &
                work(i0cmamdq+me*nsiz),work(i0cmapdq+me*nsiz), &
                work(i0cpamdq+me*nsiz),work(i0cpapdq+me*nsiz), &
                work(i0cmamdq2+me*nsiz),work(i0cmapdq2+me*nsiz), &
                work(i0cpamdq2+me*nsiz),work(i0cpapdq2+me*nsiz), &
                work(i0bmcqxxp+me*nsiz),work(i0bpcqxxp+me*nsiz), &
                work(i0bmcqxxm+me*nsiz),work(i0bpcqxxm+me*nsiz), &
                work(i0cmcqxxp+me*nsiz),work(i0cpcqxxp+me*nsiz), &
                work(i0cmcqxxm+me*nsiz),work(i0cpcqxxm+me*nsiz), &
                work(i0bmcmamdq+me*nsiz),work(i0bmcmapdq+me*nsiz), &
                work(i0bpcmamdq+me*nsiz),work(i0bpcmapdq+me*nsiz), &
                work(i0bmcpamdq+me*nsiz),work(i0bmcpapdq+me*nsiz), &
                work(i0bpcpamdq+me*nsiz),work(i0bpcpapdq+me*nsiz), &
                use_fwave,rpn3,rpt3,rptt3)
            
           cfl = dmax1(cfl,cfl1d)
            
           ! update qnew by flux differencing.
           ! (rather than maintaining arrays f, g and h for the total fluxes,
           ! the modifications are used immediately to update qnew
           ! in order to save storage.)
           
           if(index_capa == 0)then
              ! no capa array.  Standard flux differencing:
              forall (m = 1:num_eqn, i = 1:mx)
                 qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,i,me1) &
                      - dtdx * (fadd(m,i+1,me1) - fadd(m,i,me1))
              end forall
           else
              ! with capa array
              forall (m = 1:num_eqn, i = 1:mx)
                 qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,i,me1) &
                      - dtdx * (fadd(m,i+1,me1) - fadd(m,i,me1)) &
                      / aux(index_capa,i,j,k)
              end forall
           endif
           
        end do
     end do
     !$OMP END PARALLEL
    
  else if( idir == 2 )then
    
     ! perform y sweeps
     ! ==================
     !$OMP PARALLEL PRIVATE(me,me1,m,j,ma,ia,cfl1d) REDUCTION(MAX:cfl)
     me = 0 
     !$ me = omp_get_thread_num()
     me1 = me + 1 
     cfl1d = 0.d0
     
     !$OMP DO COLLAPSE(2) SCHEDULE(GUIDED)
     do k = 0, mz+1
        do i = 0, mx+1
           
           forall (m = 1:num_eqn, j = 1-num_ghost:my+num_ghost)
              ! copy data along a slice into 1d array:
              q1d(m,j,me1) = qold(m,i,j,k)
           end forall
            
           if (index_capa > 0)  then
              do j = 1-num_ghost, my+num_ghost
                 dtdy1d(j,me1) = dtdy / aux(index_capa,i,j,k)
              end do
           endif
            
           ! Since dimensional splitting is used, only aux2 is needed.
            
           if (num_aux > 0)  then
              ! There's a decent chance aux2 will all fit into cache,
              ! so keep accesses to aux as contiguous as possible.
              forall (ma=1:num_aux,ia= -1:1, j = 1-num_ghost:my+num_ghost)
                 aux2(ma, j, 2+ia,me1) = aux(ma, i+ia, j, k)
              end forall
           endif
            
           ! Store the value of i and k along this slice in the common block
           ! comxyzt in case it is needed in the Riemann solver (for
           ! variable coefficient problems)
            
           icom = i
           kcom = k
            
           ! compute modifications qadd and fadd along this slice
            
           call flux3(2,maxm,num_eqn,num_waves,num_ghost,my, &
                q1d(1,1-num_ghost,me1),dtdy1d(1-num_ghost,me1),dtdz,dtdx,dummy1,&
                aux2(1,1-num_ghost,1,me1),dummy3,num_aux, &
                method,mthlim,qadd(1,1-num_ghost,me1),fadd(1,1-num_ghost,me1),&
                gadd(1,1,-1,1-num_ghost,me1),hadd(1,1,-1,1-num_ghost,me1),cfl1d, &
                work(i0wave+me*nsiz_w),work(i0s+me*nsiz_s),work(i0amdq+me*nsiz), &
                work(i0apdq+me*nsiz),work(i0cqxx+me*nsiz), &
                work(i0bmamdq+me*nsiz),work(i0bmapdq+me*nsiz), &
                work(i0bpamdq+me*nsiz),work(i0bpapdq+me*nsiz), &
                work(i0cmamdq+me*nsiz),work(i0cmapdq+me*nsiz), &
                work(i0cpamdq+me*nsiz),work(i0cpapdq+me*nsiz), &
                work(i0cmamdq2+me*nsiz),work(i0cmapdq2+me*nsiz), &
                work(i0cpamdq2+me*nsiz),work(i0cpapdq2+me*nsiz), &
                work(i0bmcqxxp+me*nsiz),work(i0bpcqxxp+me*nsiz), &
                work(i0bmcqxxm+me*nsiz),work(i0bpcqxxm+me*nsiz), &
                work(i0cmcqxxp+me*nsiz),work(i0cpcqxxp+me*nsiz), &
                work(i0cmcqxxm+me*nsiz),work(i0cpcqxxm+me*nsiz), &
                work(i0bmcmamdq+me*nsiz),work(i0bmcmapdq+me*nsiz), &
                work(i0bpcmamdq+me*nsiz),work(i0bpcmapdq+me*nsiz), &
                work(i0bmcpamdq+me*nsiz),work(i0bmcpapdq+me*nsiz), &
                work(i0bpcpamdq+me*nsiz),work(i0bpcpapdq+me*nsiz), &
                use_fwave,rpn3,rpt3,rptt3)
            
           cfl = dmax1(cfl,cfl1d)
            
           ! update qnew by flux differencing.
           ! Note that the roles of the flux updates are changed.
           ! fadd - modifies the g-fluxes
           
           if( index_capa == 0)then
              ! no capa array.  Standard flux differencing:
              forall (m = 1:num_eqn, j = 1:my)
                 qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,j,me1) &
                      - dtdy * (fadd(m,j+1,me1) - fadd(m,j,me1))
              end forall
           else
              ! with capa array.
              forall (m = 1:num_eqn, j = 1:my)
                 qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,j,me1) &
                      - dtdy * (fadd(m,j+1,me1) - fadd(m,j,me1)) &
                      / aux(index_capa,i,j,k)
              end forall
           endif
            
        end do
     end do
     !$OMP END PARALLEL
    
  else
    
     ! perform z sweeps
     ! ==================
     !$OMP PARALLEL PRIVATE(me,me1,m,k,ma,ja,cfl1d) REDUCTION(MAX:cfl)
     me = 0 
     !$ me = omp_get_thread_num()
     me1 = me + 1 
     cfl1d = 0.d0 
     
     !$OMP DO COLLAPSE(2) SCHEDULE(GUIDED)     
     do j = 0, my+1
        do i = 0, mx+1
            
           forall (m = 1:num_eqn, k = 1-num_ghost:mz+num_ghost)
              ! copy data along a slice into 1d array:
              q1d(m,k,me1) = qold(m,i,j,k)
           end forall
            
           if (index_capa > 0)  then
              do k = 1-num_ghost, mz+num_ghost
                 dtdz1d(k,me1) = dtdz / aux(index_capa,i,j,k)
              end do
           endif
            
           ! Since dimensional splitting is used, only aux2 is needed.
            
           if (num_aux > 0)  then
              ! There's a decent chance aux2 will all fit into cache,
              ! so keep accesses to aux as contiguous as possible.
              forall (ma=1:num_aux,ja= -1:1, k = 1-num_ghost:mz+num_ghost)
                 aux2(ma, k, 2+ja,me1) = aux(ma, i, j+ja, k)
              end forall
           endif
            
           ! Store the value of i and j along this slice in the common block
           ! comxyzt in case it is needed in the Riemann solver (for
           ! variable coefficient problems)
            
           icom = i
           jcom = j
            
           ! compute modifications qadd and fadd along this slice
            
           call flux3(3,maxm,num_eqn,num_waves,num_ghost,mz, &
                q1d(1,1-num_ghost,me1),dtdz1d(1-num_ghost,me1),dtdx,dtdy,dummy1, &
                aux2(1,1-num_ghost,1,me1),dummy3,num_aux, &
                method,mthlim,qadd(1,1-num_ghost,me1),fadd(1,1-num_ghost,me1), &
                gadd(1,1,-1,1-num_ghost,me1),hadd(1,1,-1,1-num_ghost,me1),cfl1d, &
                work(i0wave+me*nsiz_w),work(i0s+me*nsiz_s),work(i0amdq+me*nsiz), &
                work(i0apdq+me*nsiz),work(i0cqxx+me*nsiz), &
                work(i0bmamdq+me*nsiz),work(i0bmapdq+me*nsiz), &
                work(i0bpamdq+me*nsiz),work(i0bpapdq+me*nsiz), &
                work(i0cmamdq+me*nsiz),work(i0cmapdq+me*nsiz), &
                work(i0cpamdq+me*nsiz),work(i0cpapdq+me*nsiz), &
                work(i0cmamdq2+me*nsiz),work(i0cmapdq2+me*nsiz), &
                work(i0cpamdq2+me*nsiz),work(i0cpapdq2+me*nsiz), &
                work(i0bmcqxxp+me*nsiz),work(i0bpcqxxp+me*nsiz), &
                work(i0bmcqxxm+me*nsiz),work(i0bpcqxxm+me*nsiz), &
                work(i0cmcqxxp+me*nsiz),work(i0cpcqxxp+me*nsiz), &
                work(i0cmcqxxm+me*nsiz),work(i0cpcqxxm+me*nsiz), &
                work(i0bmcmamdq+me*nsiz),work(i0bmcmapdq+me*nsiz), &
                work(i0bpcmamdq+me*nsiz),work(i0bpcmapdq+me*nsiz), &
                work(i0bmcpamdq+me*nsiz),work(i0bmcpapdq+me*nsiz), &
                work(i0bpcpamdq+me*nsiz),work(i0bpcpapdq+me*nsiz), &
                use_fwave,rpn3,rpt3,rptt3)
            
           cfl = dmax1(cfl,cfl1d)
            
           ! update qnew by flux differencing.
           ! Note that the roles of the flux updates are changed.
           ! fadd - modifies the h-fluxes
           
           if(index_capa == 0)then
              ! no capa array. Standard flux differencing:
              forall (m = 1:num_eqn, k = 1:mz)
                 qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,k,me1) &
                      - dtdz * (fadd(m,k+1,me1) - fadd(m,k,me1))
              end forall
           else
              ! with capa array
              forall (m = 1:num_eqn, k = 1:mz)
                 qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,k,me1) &
                      - dtdz * (fadd(m,k+1,me1) - fadd(m,k,me1)) &
                      / aux(index_capa,i,j,k)
              end forall
           endif
            
        end do
     end do
     !$OMP END PARALLEL
    
  endif

  return
end subroutine step3ds


