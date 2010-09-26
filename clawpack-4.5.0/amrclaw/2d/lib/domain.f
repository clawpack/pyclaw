c
c  ----------------------------------------------------------
c
      subroutine domain (nvar,vtime,nx,ny,naux,t0)
c
      implicit double precision (a-h,o-z)
      logical    vtime

      include  "call.i"
c
c  allocate initial coarse grid domain. set node info & initialize grid
c  initial space and time step set here too
c
      mstart = nodget(dummy)
c
c code assumes in many places that lower left corner at (0,0)
c this initial code sets the domain - assumed rectangular
c if it is too large, birect will chop it up into several rectangular
c pieces
c
      rnode(cornxlo,mstart)   = xlower
      rnode(cornylo,mstart)   = ylower
      rnode(cornyhi,mstart)   = yupper
      rnode(cornxhi,mstart)   = xupper
      node(nestlevel,mstart) = 1
      node(levelptr,mstart)  = 0
      lstart(1) = mstart

      if (tol .gt. 0.d0) then
      if (((nx/2)*2 .ne. nx) .or. (ny/2)*2 .ne. ny) then 
         write(outunit,*)" must have even number of cells"
         write(*,*)      " must have even number of cells"
         stop
      endif
      endif

      node(ndilo,mstart) = 0
      node(ndjlo,mstart) = 0
      node(ndihi,mstart) = nx-1
      node(ndjhi,mstart) = ny-1

      lfine = 1
      call  birect(mstart)
      call  ginit (mstart, .true., nvar, naux, t0)
c
c  set stable initial time step using coarse grid data
c
      if (vtime) then
         mptr = lstart(1)
         dx   = hxposs(1)
         dy   = hyposs(1)
         dt   = possk(1)
         dtgrid = dt
 60           mitot = node(ndihi,mptr)-node(ndilo,mptr) + 1 + 2*nghost
              mjtot = node(ndjhi,mptr)-node(ndjlo,mptr) + 1 + 2*nghost
              locaux = node(storeaux,mptr)
c             # added cfl to call to estdt so call.i isnt needed in estdt:
              call estdt(alloc(node(store1,mptr)),mitot,mjtot,nvar,
     1                   dx,dy,dtgrid,nghost,alloc(locaux),naux,cfl)
              dt = dmin1(dt,dtgrid)
              mptr   = node(levelptr,mptr)
            if (mptr .ne. 0) go to 60
         possk(1) = dt
      endif
c
c set rest of possk array for refined timesteps
c
      iregsz(1) = nx
      jregsz(1) = ny
      iregst(1) = 0
      jregst(1) = 0
      iregend(1) = nx-1
      jregend(1) = ny-1
      do 70 level = 2, mxnest
         iregsz(level) = iregsz(level-1) * intratx(level-1)
         jregsz(level) = jregsz(level-1) * intraty(level-1)
         possk(level)  = possk(level-1)/dble(kratio(level-1))
 70   continue
c
      return
      end
