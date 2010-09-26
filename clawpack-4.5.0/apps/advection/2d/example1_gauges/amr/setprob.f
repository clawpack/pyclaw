      subroutine setprob
      implicit double precision (a-h,o-z)
      character*25 fname
      common /comrp/ ubar,vbar
c
c     # Set the velocity for scalar advection
c     # These values are passed to the Riemann solvers rpn2.f and rpt2.f
c     # in a common block
c

c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                
      read(7,*) ubar
      read(7,*) vbar

c     # Read gauge info from setgauges.data:
      call setgauges()

      return
      end
