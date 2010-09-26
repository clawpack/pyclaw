c
c -------------------------------------------------------------
c
      subroutine stepgrid(q,fm,fp,gm,gp,mitot,mjtot,mbc,dt,dtnew,dx,dy,
     &                  nvar,xlow,ylow,time,mptr,maux,aux)
c
c          
c ::::::::::::::::::: STEPGRID ::::::::::::::::::::::::::::::::::::
c take a time step on a single grid. overwrite solution array q. 
c A modified version of the clawpack routine step2 is used.
c
c return fluxes in fm,fp and gm,gp.
c patch has room for ghost cells (mbc of them) around the grid.
c everything is the enlarged size (mitot by mjtot).
c
c mbc       = number of ghost cells  (= lwidth)
c mptr      = grid number  (for debugging)
c xlow,ylow = lower left corner of enlarged grid (including ghost cells).
c dt         = incoming time step
c dx,dy      = mesh widths for this grid
c dtnew      = return suggested new time step for this grid's soln.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      implicit double precision (a-h,o-z)
      external rpn2,rpt2

      include  "call.i"
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

      parameter (msize=max1d+4)
      parameter (mwork=msize*(maxvar*maxvar + 13*maxvar + 3*maxaux +2))

      dimension q(mitot,mjtot,nvar)
      dimension fp(mitot,mjtot,nvar),gp(mitot,mjtot,nvar)
      dimension fm(mitot,mjtot,nvar),gm(mitot,mjtot,nvar)
      dimension aux(mitot,mjtot,maux)
      dimension work(mwork)

      logical    debug,  dump
      data       debug/.false./,  dump/.false./
c
c     # set tcom = time.  This is in the common block comxyt that could
c     # be included in the Riemann solver, for example, if t is explicitly
c     # needed there.

      tcom = time

      if (dump) then
         do i = 1, mitot
         do j = 1, mjtot
            write(outunit,545) i,j,(q(i,j,ivar),ivar=1,nvar)
 545        format(2i4,4e15.7)
         end do
         end do
      endif
c
      meqn   = nvar
      mx = mitot - 2*mbc
      my = mjtot - 2*mbc
      maxm = max(mx,my)       !# size for 1d scratch array
      mbig = maxm
      xlowmbc = xlow + mbc*dx
      ylowmbc = ylow + mbc*dy

c     # method(2:7) and mthlim
c     #    are set in the amr2ez file (read by amr)
c
      method(1) = 0
c
c
c     # fluxes initialized in step2
c
      mwork0 = (maxm+2*mbc)*(12*meqn + mwaves + meqn*mwaves + 2) 
c
      if (mwork .lt. mwork0) then
         write(outunit,*) 'CLAW2 ERROR... mwork must be increased to ',
     &               mwork0
         write(*      ,*) 'CLAW2 ERROR... mwork must be increased to ',
     &               mwork0
         stop
      endif
c
c     # partition work array into pieces needed for local storage in 
c     # step2 routine. Find starting index of each piece:
c
      i0faddm = 1
      i0faddp = i0faddm + (maxm+2*mbc)*meqn
      i0gaddm = i0faddp + (maxm+2*mbc)*meqn
      i0gaddp = i0gaddm + 2*(maxm+2*mbc)*meqn
      i0q1d   = i0gaddp + 2*(maxm+2*mbc)*meqn 
      i0dtdx1 = i0q1d + (maxm+2*mbc)*meqn  
      i0dtdy1 = i0dtdx1 + (maxm+2*mbc)
      i0aux1 = i0dtdy1 + (maxm+2*mbc)
      i0aux2 = i0aux1 + (maxm+2*mbc)*maux
      i0aux3 = i0aux2 + (maxm+2*mbc)*maux
c
c
      i0next = i0aux3 + (maxm+2*mbc)*maux    !# next free space
      mused  = i0next - 1                    !# space already used
      mwork1 = mwork - mused              !# remaining space (passed to step2)

c
c
      call b4step2(mx,my,mbc,mx,my,nvar,q,
     &             xlowmbc,ylowmbc,dx,dy,time,dt,maux,aux)
c
c
c     # take one step on the conservation law:
c
      call step2(mbig,mx,my,nvar,maux,
     &           mbc,mx,my,
     &              q,aux,dx,dy,dt,cflgrid,
     &              fm,fp,gm,gp,
     &              work(i0faddm),work(i0faddp),
     &              work(i0gaddm),work(i0gaddp),
     &              work(i0q1d),work(i0dtdx1),work(i0dtdy1),
     &              work(i0aux1),work(i0aux2),work(i0aux3),
     &              work(i0next),mwork1,rpn2,rpt2)
c
c
        write(outunit,*) ' Courant # of grid ',mptr, '  is  ',cflgrid
c
        cflmax = dmax1(cflmax,cflgrid)
        cfl_level = dmax1(cfl_level,cflgrid)
c
c       # update q
        dtdx = dt/dx
        dtdy = dt/dy
        do 50 m=1,nvar
        do 50 i=mbc+1,mitot-mbc
        do 50 j=mbc+1,mjtot-mbc
         if (mcapa.eq.0) then
c
c            # no capa array.  Standard flux differencing:

           q(i,j,m) = q(i,j,m) 
     &           - dtdx * (fm(i+1,j,m) - fp(i,j,m)) 
     &           - dtdy * (gm(i,j+1,m) - gp(i,j,m)) 
         else
c            # with capa array.
           q(i,j,m) = q(i,j,m) 
     &           - (dtdx * (fm(i+1,j,m) - fp(i,j,m))
     &           + dtdy * (gm(i,j+1,m) - gp(i,j,m))) / aux(i,j,mcapa)
         endif

 50      continue
c
c
      if (method(5).eq.1) then
c        # with source term:   use Godunov splitting
         call src2(mx,my,nvar,mbc,mx,my,xlowmbc,ylowmbc,dx,dy,
     &             q,maux,aux,time,dt)
         endif
c
c
c
c     # output fluxes for debugging purposes:
      if (debug) then
         write(dbugunit,*)" fluxes for grid ",mptr
         do 830 i = mbc+1, mitot-1
            do 830 j = mbc+1, mjtot-1
               write(dbugunit,831) i,j,fm(i,j,1),fp(i,j,1),
     .                             gm(i,j,1),gp(i,j,1)
               do 830 m = 2, meqn
                  write(dbugunit,832) fm(i,j,m),fp(i,j,m),
     .            gm(i,j,m),gp(i,j,m)
  831          format(2i4,4d16.6)
  832          format(8x,4d16.6)
  830    continue
      endif

c
c
c For variable time stepping, use max speed seen on this grid to 
c choose the allowable new time step dtnew.  This will later be 
c compared to values seen on other grids.
c
       if (cflgrid .gt. 0.d0) then
           dtnew = dt*cfl/cflgrid
         else
c          # velocities are all zero on this grid so there's no 
c          # time step restriction coming from this grid.
            dtnew = rinfinity
          endif

c     # give a warning if Courant number too large...
c
      if (cflgrid .gt. cflv1) then
            write(*,810) cflgrid
            write(outunit,810) cflgrid, cflv1
  810       format('*** WARNING *** Courant number  =', d12.4,
     &              '  is larger than input cfl_max = ', d12.4)
            endif
c
      return
      end


