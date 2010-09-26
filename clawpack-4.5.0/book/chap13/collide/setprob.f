      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /cqinit/ sloc,hl,ul,hr,ur
      common /comrp/ grav
c
c     # Shallow water equations
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                
c     # Graviational constant g:
      read(7,*) grav
c
c
      return
      end
