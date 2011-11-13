c
c
c
c     ==================================================================
      subroutine step3ds(maxm,maxmx,maxmy,maxmz,meqn,mwaves,mbc,mx,my,
     &                   mz,qold,qnew,aux,dx,dy,dz,dt,method,mthlim,cfl,
     &                   qadd,fadd,gadd,hadd,q1d,dtdx1d,dtdy1d,dtdz1d,
     &                   aux1,aux2,aux3,maux,work,mwork,
     &                   rpn3,rpt3,rptt3,idir)
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
c     #       (maxm + 2*mbc)*(37*meqn + 6*maux),
c     #
c     #       when also possible reductions in flux3 are included.
c     #       However, this term is small compared to the dominating
c     #       term (maxmx+2mbc)(maxmy+2mb)*(maxmz+2mbc).
c     #-----------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      external rpn3,rpt3,rptt3
      dimension qold(meqn, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &          1-mbc:maxmz+mbc)
      dimension qnew(meqn, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &          1-mbc:maxmz+mbc)
      dimension  q1d(meqn,1-mbc:maxm+mbc)
      dimension qadd(meqn,1-mbc:maxm+mbc)
      dimension fadd(meqn,1-mbc:maxm+mbc)
      dimension gadd(meqn,2,-1:1,1-mbc:maxm+mbc)
      dimension hadd(meqn,2,-1:1,1-mbc:maxm+mbc)
      dimension aux(maux, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &              1-mbc:maxmz+mbc)
      dimension aux1(maux,1-mbc:maxm+mbc,3)
      dimension aux2(maux,1-mbc:maxm+mbc,3)
      dimension aux3(maux,1-mbc:maxm+mbc,3)
      dimension dtdx1d(1-mbc:maxm+mbc)
      dimension dtdy1d(1-mbc:maxm+mbc)
      dimension dtdz1d(1-mbc:maxm+mbc)
      dimension method(7),mthlim(mwaves)
      dimension work(mwork)
      common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom
c
c     # partition work array into pieces needed for local storage in
c     # flux2 routine.  Find starting index of each piece:
c
      i0wave     = 1
      i0s        = i0wave     + (maxm+2*mbc)*meqn*mwaves
      i0amdq     = i0s        + (maxm+2*mbc)*mwaves
      i0apdq     = i0amdq     + (maxm+2*mbc)*meqn
      i0cqxx     = i0apdq     + (maxm+2*mbc)*meqn
      i0bmamdq   = i0cqxx     + (maxm+2*mbc)*meqn
      i0bmapdq   = i0bmamdq   + (maxm+2*mbc)*meqn
      i0bpamdq   = i0bmapdq   + (maxm+2*mbc)*meqn
      i0bpapdq   = i0bpamdq   + (maxm+2*mbc)*meqn
      i0cmamdq   = i0bpapdq   + (maxm+2*mbc)*meqn
      i0cmapdq   = i0cmamdq   + (maxm+2*mbc)*meqn
      i0cpamdq   = i0cmapdq   + (maxm+2*mbc)*meqn
      i0cpapdq   = i0cpamdq   + (maxm+2*mbc)*meqn
      i0cmamdq2  = i0cpapdq   + (maxm+2*mbc)*meqn
      i0cmapdq2  = i0cmamdq2  + (maxm+2*mbc)*meqn
      i0cpamdq2  = i0cmapdq2  + (maxm+2*mbc)*meqn
      i0cpapdq2  = i0cpamdq2  + (maxm+2*mbc)*meqn
      i0bmcqxxp  = i0cpapdq2  + (maxm+2*mbc)*meqn
      i0bmcqxxm  = i0bmcqxxp  + (maxm+2*mbc)*meqn
      i0bpcqxxp  = i0bmcqxxm   + (maxm+2*mbc)*meqn
      i0bpcqxxm  = i0bpcqxxp   + (maxm+2*mbc)*meqn
      i0cmcqxxp  = i0bpcqxxm   + (maxm+2*mbc)*meqn
      i0cmcqxxm  = i0cmcqxxp   + (maxm+2*mbc)*meqn
      i0cpcqxxp  = i0cmcqxxm   + (maxm+2*mbc)*meqn
      i0cpcqxxm  = i0cpcqxxp   + (maxm+2*mbc)*meqn
      i0bmcmamdq = i0cpcqxxm   + (maxm+2*mbc)*meqn
      i0bmcmapdq = i0bmcmamdq + (maxm+2*mbc)*meqn
      i0bpcmamdq = i0bmcmapdq + (maxm+2*mbc)*meqn
      i0bpcmapdq = i0bpcmamdq + (maxm+2*mbc)*meqn
      i0bmcpamdq = i0bpcmapdq + (maxm+2*mbc)*meqn
      i0bmcpapdq = i0bmcpamdq + (maxm+2*mbc)*meqn
      i0bpcpamdq = i0bmcpapdq + (maxm+2*mbc)*meqn
      i0bpcpapdq = i0bpcpamdq + (maxm+2*mbc)*meqn
      iused      = i0bpcpapdq + (maxm+2*mbc)*meqn - 1
c
      if (iused.gt.mwork) then
c        # This shouldn't happen due to checks in claw2
         write(6,*) '*** not enough work space in step2'
         write(6,*) '*** iused = ', iused, '   mwork =',mwork
         stop
      endif
c
      mcapa = method(6)
c      maux = method(7)
      cfl = 0.d0
      dtdx = dt/dx
      dtdy = dt/dy
      dtdz = dt/dz
c
      if (mcapa.eq.0) then
c        # no capa array:
         do 5 i=1-mbc,maxm+mbc
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
            do 20 m=1,meqn
               do 20 i = 1-mbc, mx+mbc
c                 # copy data along a slice into 1d array:
                  q1d(m,i) = qold(m,i,j,k)
   20          continue
c
         if (mcapa.gt.0)  then
           do 23 i = 1-mbc, mx+mbc
               dtdx1d(i) = dtdx / aux(mcapa,i,j,k)
   23      continue
         endif
c
c        # Since dimensional splitting is used, only aux2 is needed.
c
         if (maux .gt. 0)  then
             do 22 ma=1,maux
               do 22 i = 1-mbc, mx+mbc
                 aux2(ma,i,1) = aux(ma,i,j,k-1)
                 aux2(ma,i,2) = aux(ma,i,j,k)
                 aux2(ma,i,3) = aux(ma,i,j,k+1)
   22          continue
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
            call flux3(1,maxm,meqn,mwaves,mbc,mx,
     &                 q1d,dtdx1d,dtdy,dtdz,dummy1,aux2,dummy3,maux,
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
            if(mcapa. eq. 0)then
c              # no capa array.  Standard flux differencing:
               do 30 m=1,meqn
                  do 30 i=1,mx
                     qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,i)
     &                             - dtdx * (fadd(m,i+1) - fadd(m,i))
   30             continue
            else
c              # with capa array
               do 40 m=1,meqn
                  do 40 i=1,mx
                     qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,i)
     &                        - dtdx * (fadd(m,i+1) - fadd(m,i))
     &                        / aux(mcapa,i,j,k)
   40             continue
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
            do 70 m=1,meqn
               do 70 j = 1-mbc, my+mbc
c                 # copy data along a slice into 1d array:
                  q1d(m,j) = qold(m,i,j,k)
   70          continue
c
         if (mcapa.gt.0)  then
           do 71 j = 1-mbc, my+mbc
               dtdy1d(j) = dtdy / aux(mcapa,i,j,k)
   71      continue
         endif
c
c        # Since dimensional splitting is used, only aux2 is needed.
c
         if (maux .gt. 0)  then
             do 72 ma=1,maux
               do 72 j = 1-mbc, my+mbc
                 aux2(ma,j,1) = aux(ma,i-1,j,k)
                 aux2(ma,j,2) = aux(ma,i,j,k)
                 aux2(ma,j,3) = aux(ma,i+1,j,k)
   72          continue
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
            call flux3(2,maxm,meqn,mwaves,mbc,my,
     &                 q1d,dtdy1d,dtdz,dtdx,dummy1,aux2,dummy3,maux,
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
            if( mcapa.eq. 0)then
c               # no capa array.  Standard flux differencing:
                do 80 m=1,meqn
                   do 80 j=1,my
                      qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,j)
     &                              - dtdy * (fadd(m,j+1) - fadd(m,j))
   80              continue
             else
c              #with capa array.
                do 85 m=1,meqn
                   do 85 j=1,my
                      qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,j)
     &                        - dtdy * (fadd(m,j+1) - fadd(m,j))
     &                        / aux(mcapa,i,j,k)
   85              continue
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
            do 110 m=1,meqn
               do 110 k = 1-mbc, mz+mbc
c                 # copy data along a slice into 1d array:
                  q1d(m,k) = qold(m,i,j,k)
 110           continue
c
         if (mcapa.gt.0)  then
           do 130 k = 1-mbc, mz+mbc
               dtdz1d(k) = dtdz / aux(mcapa,i,j,k)
 130       continue
         endif
c
c        # Since dimensional splitting is used, only aux2 is needed.
c
         if (maux .gt. 0)  then
             do 131 ma=1,maux
               do 131 k = 1-mbc, mz+mbc
                 aux2(ma,k,1) = aux(ma,i,j-1,k)
                 aux2(ma,k,2) = aux(ma,i,j,k)
                 aux2(ma,k,3) = aux(ma,i,j+1,k)
  131          continue
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
            call flux3(3,maxm,meqn,mwaves,mbc,mz,
     &                 q1d,dtdz1d,dtdx,dtdy,dummy1,aux2,dummy3,maux,
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
            if(mcapa .eq. 0)then
c              #no capa array. Standard flux differencing:
               do 120 m=1,meqn
                  do 120 k=1,mz
                     qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,k)
     &                             - dtdz * (fadd(m,k+1) - fadd(m,k))
 120              continue
            else
c              # with capa array
               do 145 m=1,meqn
                  do 145 k=1,mz
                     qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,k)
     &                             - dtdz * (fadd(m,k+1) - fadd(m,k))
     &                             / aux(mcapa,i,j,k)
 145              continue
            endif
c
 150  continue
c
      endif
c
      return
      end


