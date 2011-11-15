c
c     ==================================================================
      subroutine step3(maxm,maxmx,maxmy,maxmz,meqn,mwaves,mbc,mx,my,
     &                 mz,qold,qnew,aux,dx,dy,dz,dt,method,mthlim,cfl,
     &                 qadd,fadd,gadd,hadd,q1d,dtdx1d,dtdy1d,dtdz1d,
     &                 aux1,aux2,aux3,maux,work,mwork)
c     ==================================================================
c
c     # Take one time step, updating q.
c     # On entry, qold and qnew should be identical and give the
c     #    initial data for this step
c     # On exit, qnew returns values at the end of the time step.
c     #    qold is unchanged.
c
c     # qadd is used to return increments to q from flux3.
c     # fadd, gadd and hadd are used to return flux increments from flux3.
c     # See the flux3 documentation for more information.
c
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

cf2py intent(in,out) cfl
cf2py intent(in,out) qnew  
cf2py optional q1d, qadd, fadd, gadd, hadd, dtdx1d, dtdy1d, dtdz1d

      common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom

c
c
c     # partition work array into pieces needed for local storage in
c     # flux3 routine.  Find starting index of each piece:
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
c        # This shouldn't happen due to checks in claw3
         write(6,*) '*** not enough work space in step3'
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
c     # perform x-sweeps
c     ==================
c

      do 50 k = 0,mz+1
         do 50 j = 0,my+1

            forall (m = 1:meqn, i = 1-mbc:mx+mbc)
c                 # copy data along a slice into 1d array:
                  q1d(m,i) = qold(m,i,j,k)
            end forall
c
         if (mcapa.gt.0)  then
           do 23 i = 1-mbc, mx+mbc
               dtdx1d(i) = dtdx / aux(mcapa,i,j,k)
   23      continue
         endif
c
         if (maux .gt. 0)  then
            ! This odd construct may help improve cache locality.
            ! (The F95 standard says each statement in the FORALL
            ! must be executed for all indices before the next
            ! statement is started, so this is different semantically
            ! from putting all three indices in the same FORALL.)
            forall (ka = -1:1)
                forall (ma = 1:maux, i = 1-mbc:mx+mbc)
                    aux1(ma, i, 2+ka) = aux(ma, i, j-1, k+ka)
                    aux2(ma, i, 2+ka) = aux(ma, i,   j, k+ka)
                    aux3(ma, i, 2+ka) = aux(ma, i, j+1, k+ka)
                end forall
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
c           # compute modifications fadd, gadd and hadd to fluxes along
c           # this slice:
c
            call flux3(1,maxm,meqn,mwaves,mbc,mx,
     &                 q1d,dtdx1d,dtdy,dtdz,aux1,aux2,aux3,maux,
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
c
               forall (m = 1:meqn, i = 1:mx)
                  qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,i)
     &                        - dtdx * (fadd(m,i+1) - fadd(m,i))
     &                        - dtdy * (gadd(m,2,0,i) - gadd(m,1,0,i))
     &                        - dtdz * (hadd(m,2,0,i) - hadd(m,1,0,i))
                  qnew(m,i,j-1,k)   = qnew(m,i,j-1,k)
     &                              - dtdy * gadd(m,1,0,i)
     &                              - dtdz * ( hadd(m,2,-1,i)
     &                                      -   hadd(m,1,-1,i) )
                  qnew(m,i,j-1,k-1) = qnew(m,i,j-1,k-1)
     &                              - dtdy * gadd(m,1,-1,i)
     &                              - dtdz * hadd(m,1,-1,i)
                  qnew(m,i,j,k-1)   = qnew(m,i,j,k-1)
     &                              - dtdy * ( gadd(m,2,-1,i)
     &                                     -   gadd(m,1,-1,i) )
     &                              - dtdz * hadd(m,1,0,i)
                  qnew(m,i,j+1,k-1) = qnew(m,i,j+1,k-1)
     &                              + dtdy * gadd(m,2,-1,i)
     &                              - dtdz * hadd(m,1,1,i)
                  qnew(m,i,j+1,k)   = qnew(m,i,j+1,k)
     &                              + dtdy * gadd(m,2,0,i)
     &                              - dtdz * ( hadd(m,2,1,i)
     &                                     -   hadd(m,1,1,i) )
                  qnew(m,i,j+1,k+1) = qnew(m,i,j+1,k+1)
     &                              + dtdy * gadd(m,2,1,i)
     &                              + dtdz * hadd(m,2,1,i)
                  qnew(m,i,j,k+1)   = qnew(m,i,j,k+1)
     &                              - dtdy * ( gadd(m,2,1,i)
     &                                     -   gadd(m,1,1,i) )
     &                              + dtdz * hadd(m,2,0,i)
                  qnew(m,i,j-1,k+1) = qnew(m,i,j-1,k+1)
     &                              - dtdy * gadd(m,1,1,i)
     &                              + dtdz * hadd(m,2,-1,i)
               end forall
c
            else
c              # with capa array
               forall (m = 1:meqn, i = 1:mx)
                     qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,i)
     &                        - (dtdx * (fadd(m,i+1) - fadd(m,i))
     &                        +  dtdy * (gadd(m,2,0,i) - gadd(m,1,0,i))
     &                        +  dtdz * (hadd(m,2,0,i) - hadd(m,1,0,i)))
     &                        / aux(mcapa,i,j,k)

                     qnew(m,i,j-1,k)   = qnew(m,i,j-1,k)
     &                                 - (dtdy * gadd(m,1,0,i)
     &                                 +  dtdz * ( hadd(m,2,-1,i)
     &                                        -   hadd(m,1,-1,i) ))
     &                                 / aux(mcapa,i,j-1,k)
                     qnew(m,i,j-1,k-1) = qnew(m,i,j-1,k-1)
     &                                 - (dtdy * gadd(m,1,-1,i)
     &                                 +  dtdz * hadd(m,1,-1,i))
     &                                 / aux(mcapa,i,j-1,k-1)
                     qnew(m,i,j,k-1)   = qnew(m,i,j,k-1)
     &                                 - (dtdy * ( gadd(m,2,-1,i)
     &                                        -   gadd(m,1,-1,i) )
     &                                 +  dtdz * hadd(m,1,0,i))
     &                                 / aux(mcapa,i,j,k-1)
                     qnew(m,i,j+1,k-1) = qnew(m,i,j+1,k-1)
     &                                 + (dtdy * gadd(m,2,-1,i)
     &                                 -  dtdz * hadd(m,1,1,i))
     &                                 / aux(mcapa,i,j+1,k-1)
                     qnew(m,i,j+1,k)   = qnew(m,i,j+1,k)
     &                                 + (dtdy * gadd(m,2,0,i)
     &                                 - dtdz * ( hadd(m,2,1,i)
     &                                        -   hadd(m,1,1,i) ))
     &                                 / aux(mcapa,i,j+1,k)
                     qnew(m,i,j+1,k+1) = qnew(m,i,j+1,k+1)
     &                                 + (dtdy * gadd(m,2,1,i)
     &                                 +  dtdz * hadd(m,2,1,i))
     &                                 / aux(mcapa,i,j+1,k+1)
                     qnew(m,i,j,k+1)   = qnew(m,i,j,k+1)
     &                                 - (dtdy * ( gadd(m,2,1,i)
     &                                        -   gadd(m,1,1,i) )
     &                                 -  dtdz * hadd(m,2,0,i))
     &                                 / aux(mcapa,i,j,k+1)
                     qnew(m,i,j-1,k+1) = qnew(m,i,j-1,k+1)
     &                                 - (dtdy * gadd(m,1,1,i)
     &                                 -  dtdz * hadd(m,2,-1,i))
     &                                 / aux(mcapa,i,j-1,k+1)
               end forall
            endif
c
   50    continue

c
c
c     # perform y sweeps
c     ==================
c
c
      do 100 k = 0, mz+1
         do 100 i = 0, mx+1
c
            forall (m = 1:meqn, j = 1-mbc:my+mbc)
c                 # copy data along a slice into 1d array:
                  q1d(m,j) = qold(m,i,j,k)
            end forall
c
         if (mcapa.gt.0)  then
           do 71 j = 1-mbc, my+mbc
               dtdy1d(j) = dtdy / aux(mcapa,i,j,k)
   71      continue
         endif
c
         if (maux .gt. 0)  then
            ! aux1, aux2, aux3 probably fit in cache, so optimize
            ! access to aux
            forall (ma = 1:maux, ia = -1:1, j = 1-mbc:my+mbc)
                aux1(ma, j, 2+ia) = aux(ma, i+ia, j, k-1)
                aux2(ma, j, 2+ia) = aux(ma, i+ia, j,   k)
                aux3(ma, j, 2+ia) = aux(ma, i+ia, j, k+1)
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
c           # compute modifications fadd, gadd and hadd to fluxes along this
c           # slice:
c
            call flux3(2,maxm,meqn,mwaves,mbc,my,
     &                 q1d,dtdy1d,dtdz,dtdx,aux1,aux2,aux3,maux,
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
c           # gadd - modifies the h-fluxes
c           # hadd - modifies the f-fluxes
c
            if( mcapa.eq. 0)then
c               # no capa array.  Standard flux differencing:
               forall (m = 1:meqn, j = 1:my)
                  qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,j)
     &                        - dtdy * (fadd(m,j+1) - fadd(m,j))
     &                        - dtdz * (gadd(m,2,0,j) - gadd(m,1,0,j))
     &                        - dtdx * (hadd(m,2,0,j) - hadd(m,1,0,j))
                  qnew(m,i,j,k+1)   = qnew(m,i,j,k+1)
     &                              + dtdz * gadd(m,2,0,j)
     &                              - dtdx * ( hadd(m,2,1,j)
     &                                     -   hadd(m,1,1,j) )
                  qnew(m,i+1,j,k+1) = qnew(m,i+1,j,k+1)
     &                              + dtdz * gadd(m,2,1,j)
     &                              +  dtdx * hadd(m,2,1,j)
                  qnew(m,i+1,j,k)   = qnew(m,i+1,j,k)
     &                              - dtdz * ( gadd(m,2,1,j)
     &                                     -   gadd(m,1,1,j) )
     &                              + dtdx * hadd(m,2,0,j)
                  qnew(m,i+1,j,k-1) = qnew(m,i+1,j,k-1)
     &                              - dtdz * gadd(m,1,1,j)
     &                              + dtdx * hadd(m,2,-1,j)
                  qnew(m,i,j,k-1)   = qnew(m,i,j,k-1)
     &                              - dtdz * gadd(m,1,0,j)
     &                              - dtdx * ( hadd(m,2,-1,j)
     &                                     -   hadd(m,1,-1,j) )
                  qnew(m,i-1,j,k-1) = qnew(m,i-1,j,k-1)
     &                              - dtdz * gadd(m,1,-1,j)
     &                              - dtdx * hadd(m,1,-1,j)
                  qnew(m,i-1,j,k)   = qnew(m,i-1,j,k)
     &                              - dtdz * ( gadd(m,2,-1,j)
     &                                     -   gadd(m,1,-1,j) )
     &                              - dtdx * hadd(m,1,0,j)
                  qnew(m,i-1,j,k+1) = qnew(m,i-1,j,k+1)
     &                              + dtdz * gadd(m,2,-1,j)
     &                              -  dtdx*hadd(m,1,1,j)
               end forall
             else
c
c              #with capa array.
c
                forall (m = 1:meqn, j = 1:my)
                      qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,j)
     &                        - (dtdy * (fadd(m,j+1) - fadd(m,j))
     &                        +  dtdz * (gadd(m,2,0,j) - gadd(m,1,0,j))
     &                        +  dtdx * (hadd(m,2,0,j) - hadd(m,1,0,j)))
     &                        / aux(mcapa,i,j,k)
                      qnew(m,i,j,k+1)   = qnew(m,i,j,k+1)
     &                                 + (dtdz * gadd(m,2,0,j)
     &                                 -  dtdx * ( hadd(m,2,1,j)
     &                                        -   hadd(m,1,1,j) ))
     &                                 / aux(mcapa,i,j,k+1)
                      qnew(m,i+1,j,k+1) = qnew(m,i+1,j,k+1)
     &                                  + (dtdz * gadd(m,2,1,j)
     &                                  +  dtdx * hadd(m,2,1,j))
     &                                  / aux(mcapa,i+1,j,k+1)
                      qnew(m,i+1,j,k)   = qnew(m,i+1,j,k)
     &                                  - (dtdz * ( gadd(m,2,1,j)
     &                                          -   gadd(m,1,1,j) )
     &                                  -  dtdx * hadd(m,2,0,j) )
     &                                  / aux(mcapa,i+1,j,k)
                      qnew(m,i+1,j,k-1) = qnew(m,i+1,j,k-1)
     &                                  - (dtdz * gadd(m,1,1,j)
     &                                  -  dtdx * hadd(m,2,-1,j))
     &                                  / aux(mcapa,i+1,j,k-1)
                      qnew(m,i,j,k-1)   = qnew(m,i,j,k-1)
     &                                  - (dtdz * gadd(m,1,0,j)
     &                                  +  dtdx * ( hadd(m,2,-1,j)
     &                                         -   hadd(m,1,-1,j) ))
     &                                  / aux(mcapa,i,j,k-1)
                      qnew(m,i-1,j,k-1) = qnew(m,i-1,j,k-1)
     &                                  - (dtdz * gadd(m,1,-1,j)
     &                                  +  dtdx * hadd(m,1,-1,j))
     &                                  / aux(mcapa,i-1,j,k-1)
                      qnew(m,i-1,j,k)   = qnew(m,i-1,j,k)
     &                                  - (dtdz * ( gadd(m,2,-1,j)
     &                                         -   gadd(m,1,-1,j) )
     &                                  +  dtdx * hadd(m,1,0,j))
     &                                  / aux(mcapa,i-1,j,k)
                      qnew(m,i-1,j,k+1) = qnew(m,i-1,j,k+1)
     &                                  + (dtdz * gadd(m,2,-1,j)
     &                                  -  dtdx*hadd(m,1,1,j))
     &                                  / aux(mcapa,i-1,j,k+1)
                end forall
            endif
c
  100    continue
c
c
c
c     # perform z sweeps
c     ==================
c
c
      do 150 j = 0, my+1
         do 150 i = 0, mx+1
c
            forall (m = 1:meqn, k = 1-mbc:mz+mbc)
c                 # copy data along a slice into 1d array:
                  q1d(m,k) = qold(m,i,j,k)
            end forall
c
         if (mcapa.gt.0)  then
           do 130 k = 1-mbc, mz+mbc
               dtdz1d(k) = dtdz / aux(mcapa,i,j,k)
 130       continue
         endif
c
         if (maux .gt. 0)  then
            ! See the comment on the X sweeps.  This is semantically
            ! slightly different than putting all the indices in the
            ! same forall, and hopefully better for optimizing access
            ! to aux.
            forall (ja = -1:1, k = 1-mbc:mz+mbc)
                forall (ma = 1:maux)
                    aux1(ma, k, 2+ja) = aux(ma, i-1, j+ja, k)
                    aux2(ma, k, 2+ja) = aux(ma,   i, j+ja, k)
                    aux3(ma, k, 2+ja) = aux(ma, i+1, j+ja, k)
                end forall
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
c           # compute modifications fadd, gadd and hadd to fluxes along this
c           # slice:
c
            call flux3(3,maxm,meqn,mwaves,mbc,mz,
     &                 q1d,dtdz1d,dtdx,dtdy,aux1,aux2,aux3,maux,
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
c           # gadd - modifies the f-fluxes
c           # hadd - modifies the g-fluxes
c
            if(mcapa .eq. 0)then
c
c              #no capa array. Standard flux differencing:
c
               forall (m = 1:meqn, k = 1:mz)
                  qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,k)
     &                        - dtdz * (fadd(m,k+1) - fadd(m,k))
     &                        - dtdx * (gadd(m,2,0,k) - gadd(m,1,0,k))
     &                        - dtdy * (hadd(m,2,0,k) - hadd(m,1,0,k))
                  qnew(m,i,j+1,k)   = qnew(m,i,j+1,k)
     &                              - dtdx * ( gadd(m,2,1,k)
     &                                     -   gadd(m,1,1,k) )
     &                              + dtdy * hadd(m,2,0,k)
                  qnew(m,i+1,j+1,k) = qnew(m,i+1,j+1,k)
     &                              + dtdx * gadd(m,2,1,k)
     &                              + dtdy * hadd(m,2,1,k)
                  qnew(m,i+1,j,k)   = qnew(m,i+1,j,k)
     &                              + dtdx * gadd(m,2,0,k)
     &                              - dtdy * ( hadd(m,2,1,k)
     &                                     -   hadd(m,1,1,k) )
                  qnew(m,i+1,j-1,k) = qnew(m,i+1,j-1,k)
     &                              + dtdx * gadd(m,2,-1,k)
     &                              - dtdy * hadd(m,1,1,k)
                  qnew(m,i,j-1,k)   = qnew(m,i,j-1,k)
     &                              - dtdx * ( gadd(m,2,-1,k)
     &                                     -   gadd(m,1,-1,k) )
     &                              - dtdy * hadd(m,1,0,k)
                  qnew(m,i-1,j-1,k) = qnew(m,i-1,j-1,k)
     &                              - dtdx * gadd(m,1,-1,k)
     &                              - dtdy * hadd(m,1,-1,k)
                  qnew(m,i-1,j,k)   = qnew(m,i-1,j,k)
     &                              - dtdx * gadd(m,1,0,k)
     &                              - dtdy * ( hadd(m,2,-1,k)
     &                                     -   hadd(m,1,-1,k) )
                  qnew(m,i-1,j+1,k) = qnew(m,i-1,j+1,k)
     &                              - dtdx * gadd(m,1,1,k)
     &                              + dtdy * hadd(m,2,-1,k)
               end forall
            else
c
c              # with capa array
c
               forall (m = 1:meqn, k = 1:mz)
                     qnew(m,i,j,k) = qnew(m,i,j,k) + qadd(m,k)
     &                        - (dtdz * (fadd(m,k+1) - fadd(m,k))
     &                        +  dtdx * (gadd(m,2,0,k) - gadd(m,1,0,k))
     &                        +  dtdy * (hadd(m,2,0,k) - hadd(m,1,0,k)))
     &                        / aux(mcapa,i,j,k)
                     qnew(m,i,j+1,k) = qnew(m,i,j+1,k)
     &                               - (dtdx * ( gadd(m,2,1,k)
     &                                      -   gadd(m,1,1,k) )
     &                               -  dtdy * hadd(m,2,0,k))
     &                               / aux(mcapa,i,j+1,k)
                     qnew(m,i+1,j+1,k) = qnew(m,i+1,j+1,k)
     &                                 + (dtdx * gadd(m,2,1,k)
     &                                 +  dtdy * hadd(m,2,1,k))
     &                                 / aux(mcapa,i+1,j+1,k)
                     qnew(m,i+1,j,k)   = qnew(m,i+1,j,k)
     &                                 + (dtdx * gadd(m,2,0,k)
     &                                 - dtdy * ( hadd(m,2,1,k)
     &                                        -   hadd(m,1,1,k) ))
     &                                 / aux(mcapa,i+1,j,k)
                     qnew(m,i+1,j-1,k) = qnew(m,i+1,j-1,k)
     &                                 + (dtdx * gadd(m,2,-1,k)
     &                                 -  dtdy * hadd(m,1,1,k))
     &                                 / aux(mcapa,i+1,j-1,k)
                     qnew(m,i,j-1,k)   = qnew(m,i,j-1,k)
     &                                 - (dtdx * ( gadd(m,2,-1,k)
     &                                        -   gadd(m,1,-1,k) )
     &                                 +  dtdy * hadd(m,1,0,k))
     &                                 / aux(mcapa,i,j-1,k)
                     qnew(m,i-1,j-1,k) = qnew(m,i-1,j-1,k)
     &                                 - (dtdx * gadd(m,1,-1,k)
     &                                 +  dtdy * hadd(m,1,-1,k))
     &                                 / aux(mcapa,i-1,j-1,k)
                     qnew(m,i-1,j,k)   = qnew(m,i-1,j,k)
     &                                 - (dtdx * gadd(m,1,0,k)
     &                                 + dtdy * ( hadd(m,2,-1,k)
     &                                       -   hadd(m,1,-1,k) ))
     &                                 / aux(mcapa,i-1,j,k)
                     qnew(m,i-1,j+1,k) = qnew(m,i-1,j+1,k)
     &                                 - (dtdx * gadd(m,1,1,k)
     &                                 -  dtdy * hadd(m,2,-1,k))
     &                                 / aux(mcapa,i-1,j+1,k)
               end forall
         endif
c
 150  continue
c
c
      return
      end


