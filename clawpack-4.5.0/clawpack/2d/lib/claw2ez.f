c
c
c
c     =================================================================
      subroutine claw2ez(maxmx,maxmy,meqn,mwaves,mbc,maux,mwork,mthlim,
     &                   q,work,aux)
c     =================================================================
c
c     An easy-to-use clawpack driver routine for simple applications
c     Documentation is available at
c                 http://www.amath.washington.edu/~claw/doc.html
c
c     Author: Randall J. LeVeque
c     Version of August, 2002 --  CLAWPACK Version 4.1
c
      implicit double precision (a-h,o-z)
      external bc2,rpn2,rpt2,src2,b4step2

      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension work(mwork)
      dimension mthlim(mwaves)
c
      dimension method(7),dtv(5),cflv(4),nv(2),mthbc(4)
      dimension tout(100)
      logical rest
      character*12 fname
c
      common /restrt_block/ tinitial, iframe
c
c     ## New in 4.4:   Input file name changed from claw1ez.data
c     ## Open file and skip over leading lines with # comments:
      fname = 'claw.data'
      call opendatafile(55,fname)
c
      open(10,file='fort.info',status='unknown',form='formatted')
c
c
c     # Read the input in standard form from claw.data:
c
c     Number of space dimensions:  ## New in 4.4
      read(55,*) ndim
c
c     ## The remainder is unchanged from 4.3:



c     domain variables
      read(55,*) mx
      read(55,*) my
      write(6,*) '+++ mx,my: ', mx,my

c     i/o variables
      read(55,*) nout
      read(55,*) outstyle
      if (outstyle.eq.1) then
          read(55,*) tfinal
          nstepout = 1
        elseif (outstyle.eq.2) then
          read(55,*) (tout(i), i=1,nout)
          nstepout = 1
        elseif (outstyle.eq.3) then
          read(55,*) nstepout, nstop
          nout = nstop
        endif


c     timestepping variables
      read(55,*) dtv(1)
      read(55,*) dtv(2)
      read(55,*) cflv(1)
      read(55,*) cflv(2)
      read(55,*) nv(1)
c


c     # input parameters for clawpack routines
      read(55,*) method(1)
      read(55,*) method(2)
      read(55,*) method(3)
      read(55,*) method(4)
      read(55,*) method(5)
      read(55,*) method(6)
      read(55,*) method(7)

      read(55,*) meqn1
      read(55,*) mwaves1
      read(55,*) (mthlim(mw), mw=1,mwaves1)

      read(55,*) t0
      read(55,*) xlower
      read(55,*) xupper
      read(55,*) ylower
      read(55,*) yupper
c
      read(55,*) mbc1
      read(55,*) mthbc(1)
      read(55,*) mthbc(2)
      read(55,*) mthbc(3)
      read(55,*) mthbc(4)

c     # check to see if we are restarting:
      rest = .false.
c     # The next two lines may not exist in old versions of claw2ez.data.
c     # Jump over the second read statement if the 1st finds an EOF:
      read(55,*,end=199,err=199) rest
      read(55,*) iframe   !# restart from data in fort.qN file, N=iframe
 199  continue


      if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or.
     &    (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions'
         write(6,*) 'require mthbc(1) and mthbc(2) BOTH be set to 2'
         stop
         endif

      if ((mthbc(3).eq.2 .and. mthbc(4).ne.2) .or.
     &    (mthbc(4).eq.2 .and. mthbc(3).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions'
         write(6,*) 'require mthbc(3) and mthbc(4) BOTH be set to 2'
         stop
         endif

c     # These values were passed in, but check for consistency:
c
      if (method(7) .ne. maux) then
         write(6,*) '*** ERROR ***  method(7) should equal maux'
         stop
         endif
      if (meqn1 .ne. meqn) then
         write(6,*) '*** ERROR ***  meqn set wrong in input or driver'
         stop
         endif
      if (mwaves1 .ne. mwaves) then
         write(6,*) '*** ERROR ***  mwaves set wrong in input or driver'
         stop
         endif
      if (mbc1 .ne. mbc) then
         write(6,*) '*** ERROR ***  mbc set wrong in input or driver'
         stop
         endif
c
c     # check that enough storage has been allocated:
c
      if (method(5).lt.2) then
          narray = 1   !# only need one qwork array
        else
          narray = 2   !# need two qwork arrays for Strang splitting
        endif

      maxm = max0(maxmx, maxmy)
      mwork1 = (maxm+2*mbc)*(10*meqn + mwaves + meqn*mwaves
     &                      + 3*maux + 2)
     &          + narray * (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn
c
      if (mx.gt.maxmx .or. my.gt.maxmy .or. mwork.lt.mwork1) then
c        # insufficient storage
         maxmx1 = max0(mx,maxmx)
         maxmy1 = max0(my,maxmy)
         maxm1 = max0(maxmx1,maxmy1)

         mwork1 = (maxm1+2*mbc)*(10*meqn + mwaves + meqn*mwaves
     &                      + 3*maux + 2)
     &          + narray * (maxmx1 + 2*mbc) * (maxmy1 + 2*mbc) * meqn

         write(6,*) ' '
         write(6,*) '*** ERROR *** Insufficient storage allocated'
         write(6,*) 'Recompile after increasing values in driver.f:'
         write(6,611) maxmx1
         write(6,612) maxmy1
         write(6,613) mwork1
 611     format(/,'parameter (maxmx = ',i5,')')
 612     format('parameter (maxmy = ',i5,')')
 613     format('parameter (mwork = ',i9,')',/)
         stop
         endif

c
c
      write(6,*) 'running...'
      write(6,*) ' '
c
c     # grid spacing
      dx = (xupper - xlower) / float(mx)
      dy = (yupper - ylower) / float(my)
c


c     # time increments between outputing solution:
      if (outstyle .eq. 1) then
         dtout = (tfinal - t0)/float(nout)
         endif
c
c
c     # call user's routine setprob to set any specific parameters
c     # or other initialization required.
c
      call setprob
c
c     # set aux array:
c
      if (maux .gt. 0)  then
         call setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &               maux,aux)
         endif
c
c     # set initial conditions:
c
      if (rest) then
          call restart(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &          dx,dy,q)
          t0 = tinitial
        else
          call qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &          dx,dy,q,maux,aux)
          iframe = 0
        endif
c
c
      if (.not. rest) then
c        # output initial data
         call out2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,
     &          q,t0,iframe,aux,maux)
         write(6,601) iframe, t0
         endif

c
c     ----------
c     Main loop:
c     ----------
c
      tend = t0
      n0   = iframe*nstepout + 1
      do 100 n=n0,nout
         tstart = tend
         if (outstyle .eq. 1)  tend = tstart + dtout
         if (outstyle .eq. 2)  tend = tout(n)
         if (outstyle .eq. 3)  tend = tstart - 1.d0  !# single-step mode
c
         call claw2(maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &           q,aux,xlower,ylower,dx,dy,tstart,tend,dtv,
     &           cflv,nv,method,mthlim,mthbc,
     &           work,mwork,info,bc2,rpn2,rpt2,src2,b4step2)
c
c        # check to see if an error occured:
         if (info .ne. 0) then
            write(6,*) 'claw2ez aborting: Error return from claw2',
     &                 ' with info =',info
            go to 999
            endif
c
         dtv(1) = dtv(5)  !# use final dt as starting value on next call
c
c        # output solution at this time
c        ------------------------------
c
c        # if outstyle=1 or 2, then nstepout=1 and we output every time
c        # we reach this point, since claw1 was called for the entire time
c        # increment between outputs.
c
c        # if outstyle=3 then we only output if we have taken nstepout
c        # time steps since the last output.

c        # iframe is the frame number used to form file names in out1
         iframe = n/nstepout
         if (iframe*nstepout .eq. n) then
            call out2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,
     &             q,tend,iframe,aux,maux)
            write(6,601) iframe,tend
            write(10,1010) tend,info,dtv(3),dtv(4),dtv(5),
     &           cflv(3),cflv(4),nv(2)
            endif

c
c        # formats for writing out information about this call to claw:
c
  601    format('CLAW2EZ: Frame ',i4,
     &           ' output files done at time t =',
     &           d12.4,/)
c
 1010    format('tend =',d15.4,/,
     &       'info =',i5,/,'smallest dt =',d15.4,/,'largest dt =',
     &       d15.4,/,'last dt =',d15.4,/,'largest cfl =',
     &         d15.4,/,'last cfl =',d15.4,/,'steps taken =',i4,/)
c
  100    continue
c
  999 continue
c
      return
      end
