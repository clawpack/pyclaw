c     ==========================================================
      subroutine step2(maxm,num_eqn,num_waves,num_aux,num_ghost,mx,my,
     &               qold,qnew,aux,dx,dy,dt,method,mthlim,cfl,
     &               qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,use_fwave)
c     ==========================================================
c
c     # Take one time step, updating q.
c     # On entry, qold and qnew should be identical and give the
c     #    initial data for this step.
c     # On exit, qnew returns values at the end of the time step.
c     #    qold is unchanged.
c    
c     # qadd is used to return increments to q from flux2
c     # fadd and gadd are used to return flux increments from flux2.
c     # See the flux2 documentation for more information.
c
c
      implicit double precision (a-h,o-z)
      external rpn2,rpt2
      double precision qold(num_eqn,1-num_ghost:mx+num_ghost, 
     &                                  1-num_ghost:my+num_ghost)
      double precision qnew(num_eqn,1-num_ghost:mx+num_ghost, 
     &                                  1-num_ghost:my+num_ghost)
      double precision  q1d(num_eqn,1-num_ghost:maxm+num_ghost)
      double precision qadd(num_eqn,1-num_ghost:maxm+num_ghost)
      double precision fadd(num_eqn,1-num_ghost:maxm+num_ghost)
      double precision gadd(num_eqn, 2, 1-num_ghost:maxm+num_ghost)
      double precision aux(num_aux, 1-num_ghost:mx+num_ghost, 
     &                                  1-num_ghost:my+num_ghost)
      double precision aux1(num_aux, 1-num_ghost:maxm+num_ghost)
      double precision aux2(num_aux, 1-num_ghost:maxm+num_ghost)
      double precision aux3(num_aux, 1-num_ghost:maxm+num_ghost)

      double precision dtdx1d(1-num_ghost:maxm+num_ghost)
      double precision dtdy1d(1-num_ghost:maxm+num_ghost)
      integer method(7),mthlim(num_waves)
      double precision work(mwork)
      dimension qc(4)


cf2py intent(out) cfl
cf2py intent(in,out) qnew
cf2py optional q1d, qadd, fadd, gadd, dtdx1d, dtdy1d

      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
c
c
c     # partition work array into pieces needed for local storage in
c     # flux2 routine.  Find starting index of each piece:
c
      i0wave = 1
      i0s = i0wave + (maxm+2*num_ghost)*num_eqn*num_waves
      i0amdq = i0s + (maxm+2*num_ghost)*num_waves
      i0apdq = i0amdq + (maxm+2*num_ghost)*num_eqn
      i0cqxx = i0apdq + (maxm+2*num_ghost)*num_eqn
      i0bmadq = i0cqxx + (maxm+2*num_ghost)*num_eqn
      i0bpadq = i0bmadq + (maxm+2*num_ghost)*num_eqn
      iused = i0bpadq + (maxm+2*num_ghost)*num_eqn - 1
c
      if (iused.gt.mwork) then
c        # This shouldn't happen due to checks in claw2
         write(6,*) '*** not enough work space in step2'
         write(6,*) '*** iused = ', iused, '   mwork =',mwork
         stop 
      endif
c
c
      index_capa = method(6)
      num_aux = method(7)
      cfl = 0.d0
      dtdx = dt/dx
      dtdy = dt/dy
c
      if (index_capa.eq.0) then
c        # no capa array:
         do 5 i=1-num_ghost,maxm+num_ghost
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
    5    continue
      endif
    

c
c
c     # perform x-sweeps
c     ==================
c
      do 50 j = 0,my+1
c
c        # copy data along a slice into 1d arrays:
         do 21 m=1,num_eqn
            do 20 i = 1-num_ghost, mx+num_ghost
               q1d(m,i) = qold(m,i,j)
   20       continue
   21    continue
c
         if (index_capa.gt.0)  then
            do 22 i = 1-num_ghost, mx+num_ghost
               dtdx1d(i) = dtdx/aux(index_capa,i,j)
   22       continue
         endif
c
         if (num_aux .gt. 0)  then
            do 25 ma=1,num_aux
               do 24 i = 1-num_ghost, mx+num_ghost
                  aux1(ma,i) = aux(ma,i,j-1)
                  aux2(ma,i) = aux(ma,i,j  )
                  aux3(ma,i) = aux(ma,i,j+1)
   24          continue
   25       continue
         endif

c     # Store the value of j along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         jcom = j  
c           
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(1,maxm,num_eqn,num_waves,num_aux,num_ghost,mx,
     &            q1d,dtdx1d,aux1,aux2,aux3,method,mthlim,
     &            qadd,fadd,gadd,cfl1d,
     &              work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &              work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2,
     &              use_fwave)
         cfl = dmax1(cfl,cfl1d)

c        # update qnew by flux differencing.
c        # (rather than maintaining arrays f and g for the total fluxes,
c        # the modifications are used immediately to update qnew
c        # in order to save storage.)
c
         if (index_capa.eq.0) then
c
c            # no capa array.  Standard flux differencing:
            do 31 i=1,mx
               do 30 m=1,num_eqn
                  qnew(m,i,j) = qnew(m,i,j) + qadd(m,i)
     &                 - dtdx * (fadd(m,i+1) - fadd(m,i))
     &                       - dtdy * (gadd(m,2,i) - gadd(m,1,i))
                  qnew(m,i,j-1) = qnew(m,i,j-1) - dtdy * gadd(m,1,i)
                  qnew(m,i,j+1) = qnew(m,i,j+1) + dtdy * gadd(m,2,i)
   30          continue
   31       continue
c
         else
c
c            # with capa array.  
            do 41 i=1,mx
               call qcor(1,i,m,aux2,q1d,maxm,num_eqn,num_ghost,qc)
               do 40 m=1,num_eqn
                  qnew(m,i,j) = qnew(m,i,j) + qadd(m,i)
     &                 - (dtdx * (fadd(m,i+1) - fadd(m,i))
     &                         +  dtdy * (gadd(m,2,i) - gadd(m,1,i)))
     &                       / aux(index_capa,i,j)
                 qnew(m,i,j-1) = qnew(m,i,j-1) - dtdy * gadd(m,1,i)
     &                       / aux(index_capa,i,j-1)
                 qnew(m,i,j+1) = qnew(m,i,j+1) + dtdy * gadd(m,2,i)
     &                       / aux(index_capa,i,j+1)

c                # Add correction   
               qnew(m,i,j) = qnew(m,i,j)-dtdx*qc(m)/aux(index_capa,i,j)
   40          continue
   41       continue
         endif
   50 continue
c
c
c
c     # perform y sweeps
c     ==================
c
c
      do 100 i = 0, mx+1
c
c        # copy data along a slice into 1d arrays:
         do 71 m=1,num_eqn
            do 70 j = 1-num_ghost, my+num_ghost
               q1d(m,j) = qold(m,i,j)
   70       continue
   71    continue
c
         if (index_capa.gt.0)  then
            do 72 j = 1-num_ghost, my+num_ghost
               dtdy1d(j) = dtdy / aux(index_capa,i,j)
   72       continue
         endif
c
         if (num_aux .gt. 0)  then
            do 75 ma=1,num_aux
               do 74 j = 1-num_ghost, my+num_ghost
                  aux1(ma,j) = aux(ma,i-1,j)
                  aux2(ma,j) = aux(ma,i,  j)
                  aux3(ma,j) = aux(ma,i+1,j)
   74          continue
   75       continue
         endif
c
c
c     # Store the value of i along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         icom = i  
c           
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(2,maxm,num_eqn,num_waves,num_aux,num_ghost,my,
     &            q1d,dtdy1d,aux1,aux2,aux3,method,mthlim,
     &            qadd,fadd,gadd,cfl1d,
     &              work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &              work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2,
     &              use_fwave)
c
         cfl = dmax1(cfl,cfl1d)
c
c        # update qnew by flux differencing.
c        # Note that the roles of fadd and gadd are reversed for
c        # the y-sweeps -- fadd is the modification to g-fluxes and
c        # gadd is the modification to f-fluxes to the left and right.
c
         if (index_capa.eq.0) then
c
c            # no capa array.  Standard flux differencing:
            do 81 m=1,num_eqn
               do 80 j=1,my
                  qnew(m,i,j) = qnew(m,i,j) + (qadd(m,j)
     &                  - dtdy * (fadd(m,j+1) - fadd(m,j))
     &                  - dtdx * (gadd(m,2,j) - gadd(m,1,j)))
                  qnew(m,i-1,j) = qnew(m,i-1,j) - dtdx * gadd(m,1,j)
                  qnew(m,i+1,j) = qnew(m,i+1,j) + dtdx * gadd(m,2,j)
   80          continue
   81       continue
c
         else
c
c            # with capa array.  
            do 91 j=1,my
               call qcor(2,j,m,aux2,q1d,maxm,num_eqn,num_ghost,qc)
               do 90 m=1,num_eqn
                  qnew(m,i,j) = qnew(m,i,j) + qadd(m,j)
     &                  - (dtdy * (fadd(m,j+1) - fadd(m,j))
     &                    + dtdx * (gadd(m,2,j) - gadd(m,1,j)))
     &                       / aux(index_capa,i,j)
                  qnew(m,i-1,j) = qnew(m,i-1,j) - dtdx * gadd(m,1,j)
     &                       / aux(index_capa,i-1,j)
                  qnew(m,i+1,j) = qnew(m,i+1,j) + dtdx * gadd(m,2,j)
     &                       / aux(index_capa,i+1,j)

c                 # Add correction
                qnew(m,i,j) = qnew(m,i,j)-dtdy*qc(m)/aux(index_capa,i,j)
   90          continue
   91       continue
         endif
  100 continue


      return
      end
