c
c     =====================================================
      subroutine restart(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q)
c     =====================================================
c
c     # Set initial conditions for q.
c
      implicit double precision (a-h,o-z)
c
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      character*10 fname1, fname2
      common /restrt_block/ tinitial, iframe

      iunit = 16

c     # first create the file name and open file
      fname1 = 'fort.q'
     &     // char(ichar('0') + mod(iframe/1000,10)) 
     &     // char(ichar('0') + mod(iframe/100,10)) 
     &     // char(ichar('0') + mod(iframe/10,10)) 
     &     // char(ichar('0') + mod(iframe,10))
      open(iunit,file=fname1)

c     # Read grid parameters.
      read(iunit,*) igrid
      read(iunit,*) level
      read(iunit,*) mx_in
      read(iunit,*) my_in
      read(iunit,*) xlow_in
      read(iunit,*) ylow_in
      read(iunit,*) dx_in
      read(iunit,*) dy_in

c     # Test for compatibility of grid resolution.
      if (mx_in .ne. mx .or. my_in .ne. my) then
         stop 'rstart.f : data not compatible'
      endif

c     # Read variables in from old fort.qXXXX file.
      do j = 1,my
         do i = 1,mx
            read(iunit,*) (q(i,j,m),m=1,meqn)
         enddo
      enddo
      close(iunit)

c     # Read initial time in from fort.tXXXX file.      
      fname2 = 'fort.t'
     &     // char(ichar('0') + mod(iframe/1000,10)) 
     &     // char(ichar('0') + mod(iframe/100,10)) 
     &     // char(ichar('0') + mod(iframe/10,10)) 
     &     // char(ichar('0') + mod(iframe,10))
      open(iunit,file=fname2)
      read(iunit,*) tinitial
      close(iunit)


      write(*,*) 'Restarting from old output file ', fname1
      write(*,*) 'Simulation will be restarted at time t = ', tinitial
      write(*,*) 'Inital condition will not be output to a matlab ',
     &     'plot file'
      write(*,*) 

      return
      end
