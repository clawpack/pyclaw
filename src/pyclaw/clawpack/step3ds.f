c
c
c
c     ==================================================================
      subroutine step3ds(maxm,num_eqn,num_waves,num_ghost,mx,my,
     &                   mz,qold,qnew,aux,dx,dy,dz,dt,method,mthlim,cfl,
     &                   qadd,fadd,gadd,hadd,q1d,dtdx1d,dtdy1d,dtdz1d,
     &                   aux1,aux2,aux3,num_aux,work,mwork,idir)
c     ==================================================================
c
c     # Take one time step, updating q, to be used with
c     # dimensional splitting.
c     # On entry, qold and qnew should be identical and give the
c     #    initial data for this step.
c     # On exit, qnew returns values at the end of the time step.
c     #    qold is unchanged.
c
c     # This is a simplified version of the subroutine step3d
c     # Since only qadd and fadd is used in the flux updates,
c     # the solution updates below are considerably simplified
c     # compared to step3d.
c
c     #-----------------------------------------------------------------
c     # NOTE! Since dimensional splitting is used, it is possible
c     #       to reduce the memory requirement, i.e. the size of
c     #       the work array. It could be reduced with
c     #
c     #       (maxm + 2*num_ghost)*(37*num_eqn + 6*num_aux),
c     #
c     #       when also possible reductions in flux3 are included.
c     #       However, this term is small compared to the dominating
c     #       term (mx+2num_ghost)(my+2mb)*(mz+2num_ghost).
c     #-----------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      external rpn3,rpt3,rptt3
      dimension qold(num_eqn, 1-num_ghost:mx+num_ghost, 
     &          1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
      dimension qnew(num_eqn, 1-num_ghost:mx+num_ghost, 
     &          1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
      dimension  q1d(num_eqn,1-num_ghost:maxm+num_ghost)
      dimension qadd(num_eqn,1-num_ghost:maxm+num_ghost)
      dimension fadd(num_eqn,1-num_ghost:maxm+num_ghost)
      dimension gadd(num_eqn,2,-1:1,1-num_ghost:maxm+num_ghost)
      dimension hadd(num_eqn,2,-1:1,1-num_ghost:maxm+num_ghost)
      dimension aux(num_aux, 1-num_ghost:mx+num_ghost, 
     &              1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
      dimension aux1(num_aux,1-num_ghost:maxm+num_ghost,3)
      dimension aux2(num_aux,1-num_ghost:maxm+num_ghost,3)
      dimension aux3(num_aux,1-num_ghost:maxm+num_ghost,3)
      dimension dtdx1d(1-num_ghost:maxm+num_ghost)
      dimension dtdy1d(1-num_ghost:maxm+num_ghost)
      dimension dtdz1d(1-num_ghost:maxm+num_ghost)
      dimension method(7),mthlim(num_waves)
      dimension work(mwork)

cf2py intent(out) cfl
cf2py intent(in,out) qnew  
cf2py optional q1d, qadd, fadd, gadd, hadd, dtdx1d, dtdy1d, dtdz1d

      common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom
c
c     # partition work array into pieces needed for local storage in
c     # flux2 routine.  Find starting index of each piece:
c
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
c
      if (iused.gt.mwork) then
c        # This shouldn't happen due to checks in claw2
         write(6,*) '*** not enough work space in step2'
         write(6,*) '*** iused = ', iused, '   mwork =',mwork
         stop
      endif
c
      index_capa = method(6)
c      num_aux = method(7)
      cfl = 0.d0
      dtdx = dt/dx
      dtdy = dt/dy
      dtdz = dt/dz
c
      if (index_capa.eq.0) then
c        # no capa array:
         do 5 i=1-num_ghost,maxm+num_ghost
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
            dtdz1d(i) = dtdz
    5       continue
         endif
c
c
      if( idir .eq. 1)then
c
c     # perform x-sweeps
c     ==================
c
      do 50 k = 0,mz+1
         do 50 j = 0,my+1
c
            forall (m = 1:num_eqn, i = 1-num_ghost:mx+num_ghost)
c                 # copy data along a slice into 1d array:
                q1d(m,i) = qold(m,i,j,k)
            end forall
c
         if (index_capa.gt.0)  then
           do 23 i = 1-num_ghost, mx+num_ghost
               dtdx1d(i) = dtdx / aux(index_capa,i,j,k)
   23      continue
         endif
c
c        # Since dimensional splitting is used, only aux2 is needed.
c
         if (num_aux .gt. 0)  then
            forall (ma=1:num_aux,i= 1-num_ghost:mx+num_ghost, ka = -1:1)
                aux2(ma,i,2+ka) = aux(ma,i,j,k+ka)
            end forall
         endif
c
c           # Store the value of j and k along this slice in the common block
c           # comxyt in case it is needed in the Riemann solver (for
c           # variable coefficient problems)
c
            jcom = j
            kcom = k
c
c           # compute modifications qadd and fadd along this slice
c
            call flux3(1,maxm,num_eqn,num_waves,num_ghost,mx,
     &                 q1d,dtdx1d,dtdy,dtdz,dummy1,aux2,dummy3,num_aux,
     &                 method,mthlim,qadd,fadd,gadd,hadd,cfl1d,
     &                 work(i0wave),work(i0s),work(i0amdq),
     &                 work(i0apdq),work(i0cqxx),
     &                 work(i0bmamdq),work(i0bmapdq),
     &                 work(i0bpamdq),work(i0bpapdq),
     &                 work(i0cmamdq),work(i0cmapdq),
     &                 work(i0cpamdq),work(i0cpapdq),
     &                 work(i0cmamdq2),work(i0cmapdq2),
     &                 work(i0cpamdq2),work(i0cpapdq2),
     &                 work(i0bmcqxxp),work(i0bpcqxxp),
     &                 work(i0bmcqxxm),work(i0bpcqxxm),
     &                 work(i0cmcqxxp),work(i0cpcqxxp),
     &                 work(i0cmcqxxm),work(i0cpcqxxm),
     &                 work(i0bmcmamdq),work(i0bmcmapdq),
     &                 work(i0bpcmamdq),work(i0bpcmapdq),
     &                 work(i0bmcpamdq),work(i0bmcpapdq),
     &                 work(i0bpcpamdq),work(i0bpcpapdq),
     &                 rpn3,rpt3,rptt3)
c
            cfl = dmax1(cfl,cfl1d)
c
c           # update qnew by flux differencing.
c           # (rather than maintaining arrays f, g and h for the total fluxes,
c           # the modifications are used immediately to update qnew
c           # in order to save storage.)
c
            if(index_capa. eq. 0)then
c              # no capa array.  Standard flux differencing:
               forall (m = 1:num_eqn, i = 1:mx)
                     qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,i)
     &                             - dtdx * (fadd(m,i+1) - fadd(m,i))
               end forall
            else
c              # with capa array
               forall (m = 1:num_eqn, i = 1:mx)
                     qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,i)
     &                        - dtdx * (fadd(m,i+1) - fadd(m,i))
     &                        / aux(index_capa,i,j,k)
               end forall
            endif
c
   50    continue
c
      else if( idir .eq. 2 )then
c
c     # perform y sweeps
c     ==================
c
      do 100 k = 0, mz+1
         do 100 i = 0, mx+1
c
            forall (m = 1:num_eqn, j = 1-num_ghost:my+num_ghost)
c                 # copy data along a slice into 1d array:
                  q1d(m,j) = qold(m,i,j,k)
            end forall
c
         if (index_capa.gt.0)  then
           do 71 j = 1-num_ghost, my+num_ghost
               dtdy1d(j) = dtdy / aux(index_capa,i,j,k)
   71      continue
         endif
c
c        # Since dimensional splitting is used, only aux2 is needed.
c
         if (num_aux .gt. 0)  then
            ! There's a decent chance aux2 will all fit into cache,
            ! so keep accesses to aux as contiguous as possible.
            forall (ma=1:num_aux,ia= -1:1, j = 1-num_ghost:my+num_ghost)
                aux2(ma, j, 2+ia) = aux(ma, i+ia, j, k)
            end forall
         endif
c
c           # Store the value of i and k along this slice in the common block
c           # comxyzt in case it is needed in the Riemann solver (for
c           # variable coefficient problems)
c
            icom = i
            kcom = k
c
c           # compute modifications qadd and fadd along this slice
c
            call flux3(2,maxm,num_eqn,num_waves,num_ghost,my,
     &                 q1d,dtdy1d,dtdz,dtdx,dummy1,aux2,dummy3,num_aux,
     &                 method,mthlim,qadd,fadd,gadd,hadd,cfl1d,
     &                 work(i0wave),work(i0s),work(i0amdq),
     &                 work(i0apdq),work(i0cqxx),
     &                 work(i0bmamdq),work(i0bmapdq),
     &                 work(i0bpamdq),work(i0bpapdq),
     &                 work(i0cmamdq),work(i0cmapdq),
     &                 work(i0cpamdq),work(i0cpapdq),
     &                 work(i0cmamdq2),work(i0cmapdq2),
     &                 work(i0cpamdq2),work(i0cpapdq2),
     &                 work(i0bmcqxxp),work(i0bpcqxxp),
     &                 work(i0bmcqxxm),work(i0bpcqxxm),
     &                 work(i0cmcqxxp),work(i0cpcqxxp),
     &                 work(i0cmcqxxm),work(i0cpcqxxm),
     &                 work(i0bmcmamdq),work(i0bmcmapdq),
     &                 work(i0bpcmamdq),work(i0bpcmapdq),
     &                 work(i0bmcpamdq),work(i0bmcpapdq),
     &                 work(i0bpcpamdq),work(i0bpcpapdq),
     &                 rpn3,rpt3,rptt3)
c
            cfl = dmax1(cfl,cfl1d)
c
c           # update qnew by flux differencing.
c           # Note that the roles of the flux updates are changed.
c           # fadd - modifies the g-fluxes
c
            if( index_capa.eq. 0)then
c               # no capa array.  Standard flux differencing:
               forall (m = 1:num_eqn, j = 1:my)
                      qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,j)
     &                              - dtdy * (fadd(m,j+1) - fadd(m,j))
               end forall
             else
c              #with capa array.
                forall (m = 1:num_eqn, j = 1:my)
                      qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,j)
     &                        - dtdy * (fadd(m,j+1) - fadd(m,j))
     &                        / aux(index_capa,i,j,k)
                end forall
            endif
c
  100    continue
c
      else
c
c     # perform z sweeps
c     ==================
c
c
      do 150 j = 0, my+1
         do 150 i = 0, mx+1
c
            forall (m = 1:num_eqn, k = 1-num_ghost:mz+num_ghost)
c                 # copy data along a slice into 1d array:
                  q1d(m,k) = qold(m,i,j,k)
            end forall
c
         if (index_capa.gt.0)  then
           do 130 k = 1-num_ghost, mz+num_ghost
               dtdz1d(k) = dtdz / aux(index_capa,i,j,k)
 130       continue
         endif
c
c        # Since dimensional splitting is used, only aux2 is needed.
c
         if (num_aux .gt. 0)  then
            ! There's a decent chance aux2 will all fit into cache,
            ! so keep accesses to aux as contiguous as possible.
            forall (ma=1:num_aux,ja= -1:1, k = 1-num_ghost:mz+num_ghost)
                aux2(ma, k, 2+ja) = aux(ma, i, j+ja, k)
            end forall
           endif
c
c           # Store the value of i and j along this slice in the common block
c           # comxyzt in case it is needed in the Riemann solver (for
c           # variable coefficient problems)
c
            icom = i
            jcom = j
c
c           # compute modifications qadd and fadd along this slice
c
            call flux3(3,maxm,num_eqn,num_waves,num_ghost,mz,
     &                 q1d,dtdz1d,dtdx,dtdy,dummy1,aux2,dummy3,num_aux,
     &                 method,mthlim,qadd,fadd,gadd,hadd,cfl1d,
     &                 work(i0wave),work(i0s),work(i0amdq),
     &                 work(i0apdq),work(i0cqxx),
     &                 work(i0bmamdq),work(i0bmapdq),
     &                 work(i0bpamdq),work(i0bpapdq),
     &                 work(i0cmamdq),work(i0cmapdq),
     &                 work(i0cpamdq),work(i0cpapdq),
     &                 work(i0cmamdq2),work(i0cmapdq2),
     &                 work(i0cpamdq2),work(i0cpapdq2),
     &                 work(i0bmcqxxp),work(i0bpcqxxp),
     &                 work(i0bmcqxxm),work(i0bpcqxxm),
     &                 work(i0cmcqxxp),work(i0cpcqxxp),
     &                 work(i0cmcqxxm),work(i0cpcqxxm),
     &                 work(i0bmcmamdq),work(i0bmcmapdq),
     &                 work(i0bpcmamdq),work(i0bpcmapdq),
     &                 work(i0bmcpamdq),work(i0bmcpapdq),
     &                 work(i0bpcpamdq),work(i0bpcpapdq),
     &                 rpn3,rpt3,rptt3)
c
            cfl = dmax1(cfl,cfl1d)
c
c           # update qnew by flux differencing.
c           # Note that the roles of the flux updates are changed.
c           # fadd - modifies the h-fluxes
c
            if(index_capa .eq. 0)then
c              #no capa array. Standard flux differencing:
               forall (m = 1:num_eqn, k = 1:mz)
                     qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,k)
     &                             - dtdz * (fadd(m,k+1) - fadd(m,k))
               end forall
            else
c              # with capa array
               forall (m = 1:num_eqn, k = 1:mz)
                     qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,k)
     &                             - dtdz * (fadd(m,k+1) - fadd(m,k))
     &                             / aux(index_capa,i,j,k)
               end forall
            endif
c
 150  continue
c
      endif
c
      return
      end


