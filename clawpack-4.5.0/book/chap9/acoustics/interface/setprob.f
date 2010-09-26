      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /cqinit/ beta,ic
      common /comaux/ rhol,cl,rhor,cr
      common /comlim/ mylim,mrplim(2)
c
c     # Set the material parameters for the acoustic equations
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                
c     # choice of initial data:
      read(7,*) ic
c     # beta for initial conditions:
      read(7,*) beta
c
c     # Piecewise constant medium with single interface at x=0
c     # Density and sound speed to left and right:
      read(7,*) rhol
      read(7,*) cl
      read(7,*) rhor
      read(7,*) cr

      return
      end
