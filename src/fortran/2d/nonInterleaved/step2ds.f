c
c
c
c
c     ==========================================================
      subroutine step2ds(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &               qold,qnew,aux,dx,dy,dt,method,mthlim,cfl,
     &               qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,rpn2,rpt2,ids)
c     ==========================================================
c
c     # Take one time step, updating q.
c     # On entry, qold and qnew should be identical and give the
c     #    initial data for this step
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
      dimension qold(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension qnew(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension  q1d(1-mbc:maxm+mbc, meqn)
      dimension qadd(1-mbc:maxm+mbc, meqn)
      dimension fadd(1-mbc:maxm+mbc, meqn)
      dimension gadd(1-mbc:maxm+mbc, meqn, 2)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, *)
      dimension aux1(1-mbc:maxm+mbc, *)
      dimension aux2(1-mbc:maxm+mbc, *)
      dimension aux3(1-mbc:maxm+mbc, *)

      dimension dtdx1d(1-mbc:maxm+mbc)
      dimension dtdy1d(1-mbc:maxm+mbc)
      dimension method(7),mthlim(mwaves)
      dimension work(mwork)

      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
c
c
c     # partition work array into pieces needed for local storage in
c     # flux2 routine.  Find starting index of each piece:
c
      i0wave = 1
      i0s = i0wave + (maxm+2*mbc)*meqn*mwaves
      i0amdq = i0s + (maxm+2*mbc)*mwaves
      i0apdq = i0amdq + (maxm+2*mbc)*meqn
      i0cqxx = i0apdq + (maxm+2*mbc)*meqn
      i0bmadq = i0cqxx + (maxm+2*mbc)*meqn
      i0bpadq = i0bmadq + (maxm+2*mbc)*meqn
      iused = i0bpadq + (maxm+2*mbc)*meqn - 1
c
      if (iused.gt.mwork) then
c        # This shouldn't happen due to checks in claw2
         write(6,*) '*** not enough work space in step2'
         write(6,*) '*** iused = ', iused, '   mwork =',mwork
         stop 
      endif
c
c
      mcapa = method(6)
      maux = method(7)
      cfl = 0.d0
      dtdx = dt/dx
      dtdy = dt/dy
c
      if (mcapa.eq.0) then
c        # no capa array:
         do 5 i=1-mbc,maxm+mbc
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
    5    continue
      endif
c
      if( ids.eq.1 )then
c
c     # perform x-sweeps
c     ==================
c
c     # note that for dimensional splitting we sweep over the rows of
c     # ghosts cells as well as the interior.  This updates the ghost
c     # cell values to the intermediate state as needed in the following 
c     # sweep in the y-direction.
c
      do 50 j = 1-mbc,my+mbc
c
c        # copy data along a slice into 1d arrays:
         do 21 m=1,meqn
            do 20 i = 1-mbc, mx+mbc
               q1d(i,m) = qold(i,j,m)
   20       continue
   21    continue
c
         if (mcapa.gt.0)  then
            do 22 i = 1-mbc, mx+mbc
               dtdx1d(i) = dtdx / aux(i,j,mcapa)
   22       continue
         endif
c
         if (maux .gt. 0)  then
             do 23 ma=1,maux
               do 23 i = 1-mbc, mx+mbc
                 aux2(i,ma) = aux(i,j  ,ma)
   23          continue
c
             if(j .ne. 1-mbc)then
                do 24 ma=1,maux
                   do 24 i = 1-mbc, mx+mbc
                      aux1(i,ma) = aux(i,j-1,ma)
   24              continue
                endif
c
             if(j .ne. my+mbc)then
                do 25 ma=1,maux
                   do 25 i = 1-mbc, mx+mbc
                      aux3(i,ma) = aux(i,j+1,ma)
   25              continue
                endif
c
             endif
c
c        # Store the value of j along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         jcom = j  
c           
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(1,maxm,meqn,mwaves,mbc,mx,
     &            q1d,dtdx1d,aux1,aux2,aux3,method,mthlim,
     &            qadd,fadd,gadd,cfl1d,
     &              work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &              work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2)
         cfl = dmax1(cfl,cfl1d)
c
c        # update qnew by flux differencing.
c        # (rather than maintaining arrays f and g for the total fluxes,
c        # the modifications are used immediately to update qnew
c        # in order to save storage.)
c
         if (mcapa.eq.0) then
c
c            # no capa array.  Standard flux differencing:
            do 31 m=1,meqn
               do 30 i=1,mx
                  qnew(i,j,m) = qnew(i,j,m) + qadd(i,m)
     &                 - dtdx * (fadd(i+1,m) - fadd(i,m))
   30          continue
   31       continue
c
         else
c
c            # with capa array.  
            do 41 m=1,meqn
               do 40 i=1,mx
                  qnew(i,j,m) = qnew(i,j,m) + qadd(i,m)
     &                        - dtdx * (fadd(i+1,m) - fadd(i,m))
     &                        / aux(i,j,mcapa)
   40          continue
   41       continue
         endif
   50 continue
c
      endif
c
      if( ids.eq.2 )then
c
c     # perform y sweeps
c     ==================
c
c
      do 100 i = 1-mbc, mx+mbc
c
c        # copy data along a slice into 1d arrays:
         do 71 m=1,meqn
            do 70 j = 1-mbc, my+mbc
               q1d(j,m) = qold(i,j,m)
   70       continue
   71    continue
c
         if (mcapa.gt.0)  then
            do 72 j = 1-mbc, my+mbc
               dtdy1d(j) = dtdy / aux(i,j,mcapa)
   72       continue
         endif
c
         if (maux .gt. 0)  then
c
             do 73 ma=1,maux
               do 73 j = 1-mbc, my+mbc
                 aux2(j,ma) = aux(i,j,ma)
   73          continue
c
             if(i .ne. 1-mbc)then
                do 74 ma=1,maux
                   do 74 j = 1-mbc, my+mbc
                      aux1(j,ma) = aux(i-1,j,ma)
   74              continue
                endif
c
             if(i .ne. mx+mbc)then
                do 75 ma=1,maux
                   do 75 j = 1-mbc, my+mbc
                      aux3(j,ma) = aux(i+1,j,ma)
   75              continue
                endif
c
             endif            
c
c     # Store the value of i along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         icom = i  
c           
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(2,maxm,meqn,mwaves,mbc,my,
     &            q1d,dtdy1d,aux1,aux2,aux3,method,mthlim,
     &            qadd,fadd,gadd,cfl1d,
     &              work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &              work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2)
c
         cfl = dmax1(cfl,cfl1d)
c
c        # update qnew by flux differencing.
c        # Note that the roles of fadd and gadd are reversed for
c        # the y-sweeps -- fadd is the modification to g-fluxes and
c        # gadd is the modification to f-fluxes to the left and right.
c
         if (mcapa.eq.0) then
c
c            # no capa array.  Standard flux differencing:
            do 81 m=1,meqn
               do 80 j=1,my
                  qnew(i,j,m) = qnew(i,j,m) + qadd(j,m)
     &                  - dtdy * (fadd(j+1,m) - fadd(j,m))
   80          continue
   81       continue
c
         else
c
c            # with capa array.  
            do 91 m=1,meqn
               do 90 j=1,my
                  qnew(i,j,m) = qnew(i,j,m) + qadd(j,m)
     &                  - dtdy * (fadd(j+1,m) - fadd(j,m))
     &                        / aux(i,j,mcapa)
   90          continue
   91       continue
         endif
  100 continue
c
      endif
c
c
      return
      end
